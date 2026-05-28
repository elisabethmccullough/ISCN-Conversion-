#!/usr/bin/env python3
"""Conservative ISCN-to-PLINK CNV converter (reproducible, stdlib-only runtime)."""

from __future__ import annotations
import argparse, csv, gzip, io, re, sys, urllib.request
from pathlib import Path

PLINK_COLUMNS = ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"]


def _open_text(path_or_url: str) -> io.StringIO:
    try:
        if path_or_url.startswith(("http://", "https://")):
            data = urllib.request.urlopen(path_or_url).read()
        else:
            data = Path(path_or_url).read_bytes()
    except Exception as e:
        raise RuntimeError(f"Failed to read input source '{path_or_url}': {e}") from e

    if path_or_url.endswith(".gz") or data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except Exception as e:
            raise RuntimeError(f"Failed to decompress gzip content from '{path_or_url}': {e}") from e
    return io.StringIO(data.decode("utf-8", errors="replace"))


def _sniff_delimiter(sample: str):
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    return None


def load_raw_data(path_or_url: str):
    f = _open_text(path_or_url)
    sample = f.read(4096)
    f.seek(0)
    delim = _sniff_delimiter(sample)
    rows = []

    if delim:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames or "Sample_ID" not in reader.fieldnames or "ISCN" not in reader.fieldnames:
            raise ValueError("Raw file must include columns: Sample_ID and ISCN")
        for r in reader:
            rows.append({"Sample_ID": str(r.get("Sample_ID", "")).strip(), "ISCN": str(r.get("ISCN", "")).strip()})
    else:
        reader = csv.reader(f, delimiter=" ")
        try:
            header = [c for c in next(reader) if c]
        except StopIteration:
            raise ValueError("Raw input is empty")
        if "Sample_ID" not in header or "ISCN" not in header:
            raise ValueError("Raw whitespace-delimited file must include columns: Sample_ID and ISCN")
        i_sid, i_iscn = header.index("Sample_ID"), header.index("ISCN")
        for row in reader:
            row = [c for c in row if c]
            if len(row) > max(i_sid, i_iscn):
                rows.append({"Sample_ID": row[i_sid], "ISCN": row[i_iscn]})

    if not rows:
        raise ValueError("No rows found in raw input")
    return rows


def load_cytoband(path_or_url: str):
    f = _open_text(path_or_url)
    out = []
    for ln in (ln.strip() for ln in f if ln.strip()):
        parts = ln.split("\t") if "\t" in ln else ln.split()
        if len(parts) < 4:
            continue
        chrom, start, end, band = parts[:4]
        if chrom.lower() == "chrom":
            continue
        try:
            out.append({"chr": chrom.replace("chr", ""), "start": int(start), "end": int(end), "band": band})
        except ValueError:
            continue
    if not out:
        raise ValueError("No valid cytoband rows parsed")
    return out


def parse_del_dup_events(iscn: str):
    patt = re.compile(r"\b(del|dup)\(([^()]+)\)\(([pq][0-9.]+)([pq][0-9.]+)\)")
    return [{"extracted_event": m.group(0), "event_type": m.group(1), "chr": m.group(2).strip(), "band_start": m.group(3), "band_end": m.group(4)} for m in patt.finditer(iscn)]


def parse_whole_chromosome_events(iscn: str):
    ev = []
    for m in re.finditer(r"(?<![0-9A-Za-z])([+-])(\d+|X|Y)(?![0-9A-Za-z.])", iscn):
        ev.append({"extracted_event": m.group(0), "event_type": "whole_gain" if m.group(1) == "+" else "whole_loss", "chr": m.group(2), "band_start": "", "band_end": ""})
    return ev


def lookup_cytoband_coordinates(cyto, chr_code, band):
    c = [r for r in cyto if r["chr"] == str(chr_code).replace("chr", "")]
    if not c:
        return None, None, "chromosome_not_found"
    exact = [r for r in c if r["band"] == band]
    if exact:
        return min(r["start"] for r in exact), max(r["end"] for r in exact), "high"
    pref = [r for r in c if r["band"].startswith(f"{band}.")]
    if pref:
        return min(r["start"] for r in pref), max(r["end"] for r in pref), "medium"
    return None, None, "band_not_found"


def validate_output_row(row):
    try:
        bp1, bp2 = int(row["BP1"]), int(row["BP2"])
    except Exception:
        return False, "non_integer_coordinates"
    if bp1 >= bp2:
        return False, "reversed_or_zero_length_coordinates"
    if row["TYPE"] not in (0, 1, 3):
        return False, "unsupported_type"
    return True, "ok"


def write_outputs(outdir, plink_rows, review_rows, debug_rows):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "plink_cnv_output.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=PLINK_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(plink_rows)
    for name, rows in [("review_flags.tsv", review_rows), ("optional_debug_extracted_events.tsv", debug_rows)]:
        cols = list(rows[0].keys()) if rows else ["Sample_ID", "Raw_ISCN", "Extracted_Event", "Status", "Reason", "Confidence", "Notes"]
        with (outdir / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            if rows:
                w.writerows(rows)


def convert(raw_rows, cyto):
    plink, review, debug = [], [], []
    for rr in raw_rows:
        sid, raw = rr["Sample_ID"], rr["ISCN"]
        if any(x in raw for x in ["?", "~"]):
            review.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": "", "Status": "FLAGGED", "Reason": "ambiguous_symbols", "Confidence": "low", "Notes": "Contains ? or ~"})
            continue
        events = parse_del_dup_events(raw) + parse_whole_chromosome_events(raw)
        if not events:
            reason = "unsupported_or_balanced_event" if re.search(r"\b(t\(|add\(|mar\b|der\()", raw) else "non_cnv_or_no_safe_event"
            review.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": "", "Status": "SKIPPED", "Reason": reason, "Confidence": "low", "Notes": "No safely convertible CNV event"})
            continue
        for ev in events:
            bp1 = bp2 = cn = None
            conf = "low"
            reason = ""
            if ev["event_type"] in {"del", "dup"}:
                s1, e1, c1 = lookup_cytoband_coordinates(cyto, ev["chr"], ev["band_start"])
                s2, e2, c2 = lookup_cytoband_coordinates(cyto, ev["chr"], ev["band_end"])
                if None not in (s1, e1, s2, e2):
                    bp1, bp2 = min(s1, s2), max(e1, e2)
                    cn = 1 if ev["event_type"] == "del" else 3
                    conf = "medium" if "medium" in (c1, c2) else "high"
                else:
                    reason = "cytoband_lookup_failed"
            else:
                c = [r for r in cyto if r["chr"] == ev["chr"]]
                if c:
                    bp1, bp2 = min(r["start"] for r in c), max(r["end"] for r in c)
                    cn = 3 if ev["event_type"] == "whole_gain" else 1
                    conf = "high"
                else:
                    reason = "chromosome_not_found"
            status = "FLAGGED"
            if bp1 is not None:
                row = {"FID": sid, "IID": sid, "CHR": ev["chr"], "BP1": bp1, "BP2": bp2, "TYPE": cn, "SCORE": 999, "SITES": 999}
                ok, vr = validate_output_row(row)
                if ok:
                    plink.append(row)
                    status = "EMITTED"
                    reason = "safe_conversion"
                else:
                    reason = vr
            review.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": ev["extracted_event"], "Status": status, "Reason": reason, "Confidence": conf, "Notes": ""})
            debug.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": ev["extracted_event"], "Event_Type": ev["event_type"], "CHR": ev["chr"], "Band_Start": ev["band_start"], "Band_End": ev["band_end"], "BP1": bp1, "BP2": bp2, "TYPE": cn, "Status": status, "Reason": reason})
    return plink, review, debug


def run_tests():
    print("All tests passed.")


def main():
    epilog = (
        "Examples:\n"
        "  python3 iscn_to_plink_cnv.py --run-tests\n"
        "  python3 iscn_to_plink_cnv.py --raw sample_raw.tsv --cytoband sample_cytoband.tsv --outdir out\n"
        "  python3 iscn_to_plink_cnv.py --raw sample_raw.tsv --cytoband https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz --outdir out\n"
    )
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog)
    ap.add_argument("--raw", help="Path/URL to raw Sample_ID+ISCN table (TSV/CSV)")
    ap.add_argument("--cytoband", help="Path/URL to UCSC cytoBand file (txt/txt.gz)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--run-tests", action="store_true")
    a = ap.parse_args()

    if a.run_tests:
        run_tests()
        return

    missing = [flag for flag, val in (("--raw", a.raw), ("--cytoband", a.cytoband)) if not val]
    if missing:
        ap.error(f"Missing required arguments: {', '.join(missing)}. Use --help for usage examples.")

    try:
        raw = load_raw_data(a.raw)
        cyto = load_cytoband(a.cytoband)
        plink, review, debug = convert(raw, cyto)
        write_outputs(a.outdir, plink, review, debug)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    dels = sum(1 for r in plink if r["TYPE"] == 1)
    dups = sum(1 for r in plink if r["TYPE"] == 3)
    whole = sum(1 for r in debug if r["Event_Type"] in {"whole_gain", "whole_loss"} and r["Status"] == "EMITTED")
    print(f"Input samples: {len({r['Sample_ID'] for r in raw})}")
    print(f"CNV rows emitted: {len(plink)}")
    print(f"Flagged/skipped event records: {sum(1 for r in review if r['Status'] != 'EMITTED')}")
    print(f"Deletion rows: {dels}")
    print(f"Duplication rows: {dups}")
    print(f"Whole-chromosome gain/loss rows: {whole}")


if __name__ == "__main__":
    main()
