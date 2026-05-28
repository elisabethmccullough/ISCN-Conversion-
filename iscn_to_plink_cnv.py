#!/usr/bin/env python3
"""Conservative ISCN-to-PLINK CNV converter (reproducible, stdlib-only runtime)."""

from __future__ import annotations
import argparse, csv, gzip, io, re, sys, tempfile, urllib.request
from pathlib import Path

TYPE_COLUMN = "TYPE (# of copies)"
PLINK_COLUMNS = ["FID", "IID", "CHR", "BP1", "BP2", TYPE_COLUMN, "SCORE", "SITES"]
REVIEW_COLUMNS = ["Sample_ID", "Raw_ISCN", "Extracted_Event", "Status", "Reason", "Confidence", "Notes"]
DEBUG_COLUMNS = ["Sample_ID", "Raw_ISCN", "Extracted_Event", "Event_Type", "CHR", "Band_Start", "Band_End", "BP1", "BP2", TYPE_COLUMN, "Status", "Reason"]


def _open_text(path_or_url: str) -> io.StringIO:
    try:
        if path_or_url.startswith(("http://", "https://")):
            request = urllib.request.Request(path_or_url, headers={"User-Agent": "iscn-to-plink-cnv/1.0"})
            data = urllib.request.urlopen(request).read()
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
    # Safe prefix expansion handles both broad integer bands (q13 -> q13.1,
    # q13.2, ...) and decimal parent bands (p11.2 -> p11.21, p11.22, ...).
    prefix = f"{band}." if "." not in band else band
    pref = [r for r in c if r["band"].startswith(prefix)]
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
    if not str(row.get("CHR", "")).strip():
        return False, "missing_chromosome"
    if row.get(TYPE_COLUMN) not in (1, 3):
        return False, "unsupported_type"
    return True, "ok"




def validate_plink_rows(plink_rows):
    """Validate the complete PLINK table before writing it to disk."""
    seen = set()
    for i, row in enumerate(plink_rows, start=1):
        if list(row.keys()) != PLINK_COLUMNS:
            raise ValueError(f"PLINK row {i} columns do not match expected columns: {PLINK_COLUMNS}")
        if not str(row["FID"]).strip() or not str(row["IID"]).strip():
            raise ValueError(f"PLINK row {i} has missing FID or IID")
        if row["FID"] != row["IID"]:
            raise ValueError(f"PLINK row {i} has FID and IID mismatch")
        ok, reason = validate_output_row(row)
        if not ok:
            raise ValueError(f"PLINK row {i} failed validation: {reason}")
        key = tuple(row[col] for col in PLINK_COLUMNS)
        if key in seen:
            raise ValueError(f"PLINK row {i} duplicates an earlier row")
        seen.add(key)


def validate_written_plink_tsv(path):
    """Read the written TSV back and confirm it is a clean 8-column table.

    The preferred read-back test uses pandas with sep="\t" when pandas is
    installed. In minimal environments without pandas, the same structural
    checks are performed with the standard-library csv module so --run-tests
    remains runnable.
    """
    path = Path(path)
    forbidden_markers = ("git apply", "diff --git", "/dev/null", "EOF")
    text = path.read_text()
    if any(marker in text for marker in forbidden_markers):
        raise ValueError("plink_cnv_output.tsv contains patch/diff text instead of only table rows")

    with path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError("plink_cnv_output.tsv is empty")
    if rows[0] != PLINK_COLUMNS:
        raise ValueError(f"Unexpected TSV header: {rows[0]}")
    for i, row in enumerate(rows[1:], start=2):
        if len(row) != len(PLINK_COLUMNS):
            raise ValueError(f"TSV row {i} has {len(row)} columns; expected {len(PLINK_COLUMNS)}")
        if not row[0] or not row[1] or row[0] != row[1]:
            raise ValueError(f"TSV row {i} has invalid FID/IID values")

    try:
        import pandas as pd  # type: ignore
    except ImportError:
        return "csv"

    df = pd.read_csv(path, sep="\t")
    if list(df.columns) != PLINK_COLUMNS:
        raise ValueError(f"pandas read-back columns do not match expected columns: {list(df.columns)}")
    if not (df["FID"].notna().all() and df["IID"].notna().all()):
        raise ValueError("pandas read-back found missing FID/IID values")
    if not (df["FID"].astype(str) == df["IID"].astype(str)).all():
        raise ValueError("pandas read-back found FID/IID mismatches")
    if not (df["BP1"].astype(int) < df["BP2"].astype(int)).all():
        raise ValueError("pandas read-back found BP1 >= BP2")
    if not df[TYPE_COLUMN].isin([1, 3]).all():
        raise ValueError("pandas read-back found unsupported TYPE (# of copies)")
    if df.duplicated().any():
        raise ValueError("pandas read-back found duplicate rows")
    return "pandas"


def write_outputs(outdir, plink_rows, review_rows, debug_rows):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    validate_plink_rows(plink_rows)
    with (outdir / "plink_cnv_output.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=PLINK_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(plink_rows)

    # Convenience CSV copy for users who want to download/open the main
    # PLINK-style output directly in spreadsheet software. The TSV remains
    # the canonical required output file.
    with (outdir / "plink_cnv_output.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=PLINK_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(plink_rows)
    for name, rows in [("review_flags.tsv", review_rows), ("optional_debug_extracted_events.tsv", debug_rows)]:
        cols = DEBUG_COLUMNS if name.startswith("optional_debug") else REVIEW_COLUMNS
        with (outdir / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            if rows:
                w.writerows(rows)


def _review_row(sample_id, raw_iscn, extracted_event, status, reason, confidence="low", notes=""):
    """Build one standardized review row."""
    return {
        "Sample_ID": sample_id,
        "Raw_ISCN": raw_iscn,
        "Extracted_Event": extracted_event,
        "Status": status,
        "Reason": reason,
        "Confidence": confidence,
        "Notes": notes,
    }


def _unsupported_tokens(iscn: str):
    """Return unsupported structural-event tokens that should be reviewed, not inferred as CNVs."""
    tokens = []
    for pattern in (r"\bt\([^)]*\)(?:\([^)]*\))?", r"\bder\([^)]*\)(?:\([^)]*\))?", r"\badd\([^)]*\)(?:\([^)]*\))?", r"\bmar\b"):
        tokens.extend(m.group(0) for m in re.finditer(pattern, iscn))
    return tokens


def _make_debug_row(sample_id, raw_iscn, event, bp1, bp2, copy_number, status, reason):
    """Build one standardized parser/debug row."""
    return {
        "Sample_ID": sample_id,
        "Raw_ISCN": raw_iscn,
        "Extracted_Event": event["extracted_event"],
        "Event_Type": event["event_type"],
        "CHR": event["chr"],
        "Band_Start": event["band_start"],
        "Band_End": event["band_end"],
        "BP1": min(bp1, bp2) if bp1 is not None and bp2 is not None else bp1,
        "BP2": max(bp1, bp2) if bp1 is not None and bp2 is not None else bp2,
        TYPE_COLUMN: copy_number,
        "Status": status,
        "Reason": reason,
    }


def convert(raw_rows, cyto):
    """Convert clear CNV events while routing uncertain/unsupported content to review files."""
    plink, review, debug = [], [], []
    emitted_keys = set()

    for rr in raw_rows:
        sid, raw = rr["Sample_ID"], rr["ISCN"]

        # Question marks and approximate ranges make the whole ISCN unsafe to convert.
        if any(x in raw for x in ["?", "~"]):
            review.append(_review_row(sid, raw, "", "FLAGGED", "ambiguous_symbols", notes="Contains ? or ~"))
            continue

        events = parse_del_dup_events(raw) + parse_whole_chromosome_events(raw)
        unsupported = _unsupported_tokens(raw)
        for token in unsupported:
            review.append(_review_row(sid, raw, token, "FLAGGED", "unsupported_or_balanced_event", notes="Not converted; clear CNV events in same ISCN are still processed separately."))

        if not events:
            reason = "unsupported_or_balanced_event" if unsupported else "non_cnv_or_no_safe_event"
            review.append(_review_row(sid, raw, "", "SKIPPED", reason, notes="No safely convertible CNV event"))
            continue

        for ev in events:
            bp1 = bp2 = copy_number = None
            confidence = "low"
            reason = ""

            if ev["event_type"] in {"del", "dup"}:
                start1, end1, conf1 = lookup_cytoband_coordinates(cyto, ev["chr"], ev["band_start"])
                start2, end2, conf2 = lookup_cytoband_coordinates(cyto, ev["chr"], ev["band_end"])
                if None not in (start1, end1, start2, end2):
                    # Use the start of the first band and end of the second band, then sort.
                    bp1, bp2 = sorted((int(start1), int(end2)))
                    copy_number = 1 if ev["event_type"] == "del" else 3
                    confidence = "medium" if "medium" in (conf1, conf2) else "high"
                else:
                    reason = "cytoband_lookup_failed"
            else:
                chrom_rows = [r for r in cyto if r["chr"] == ev["chr"]]
                if chrom_rows:
                    bp1, bp2 = min(r["start"] for r in chrom_rows), max(r["end"] for r in chrom_rows)
                    copy_number = 3 if ev["event_type"] == "whole_gain" else 1
                    confidence = "high"
                else:
                    reason = "chromosome_not_found"

            status = "FLAGGED"
            if bp1 is not None and bp2 is not None and copy_number is not None:
                chr_code = str(ev["chr"]).replace("chr", "")
                bp1, bp2 = sorted((int(bp1), int(bp2)))
                row = {
                    "FID": sid,
                    "IID": sid,
                    "CHR": chr_code,
                    "BP1": bp1,
                    "BP2": bp2,
                    "TYPE (# of copies)": copy_number,
                    "SCORE": 999,
                    "SITES": 999,
                }
                ok, validation_reason = validate_output_row(row)
                dedupe_key = tuple(row[col] for col in PLINK_COLUMNS)
                if ok and dedupe_key not in emitted_keys:
                    plink.append(row)
                    emitted_keys.add(dedupe_key)
                    status = "EMITTED"
                    reason = "safe_conversion"
                elif ok:
                    status = "SKIPPED"
                    reason = "duplicate_output_row"
                else:
                    reason = validation_reason

            review.append(_review_row(sid, raw, ev["extracted_event"], status, reason, confidence))
            debug.append(_make_debug_row(sid, raw, ev, bp1, bp2, copy_number, status, reason))

    return plink, review, debug


def _fake_cytobands():
    """Small hg38-like cytoband table for built-in tests."""
    return [
        {"chr": "5", "start": 67000000, "end": 70000000, "band": "q13"},
        {"chr": "5", "start": 150000000, "end": 160500000, "band": "q33"},
        {"chr": "6", "start": 0, "end": 2300000, "band": "p25.3"},
        {"chr": "6", "start": 25000000, "end": 30500000, "band": "p22.1"},
        {"chr": "6", "start": 105000000, "end": 110000000, "band": "q21"},
        {"chr": "6", "start": 158000000, "end": 160600000, "band": "q25.3"},
        {"chr": "7", "start": 0, "end": 159000000, "band": "p11"},
        {"chr": "7", "start": 77900000, "end": 80000000, "band": "q21"},
        {"chr": "7", "start": 120000000, "end": 127500000, "band": "q31.33"},
        {"chr": "8", "start": 0, "end": 145138636, "band": "p11"},
        {"chr": "10", "start": 0, "end": 133797422, "band": "p11"},
        {"chr": "12", "start": 0, "end": 3200000, "band": "p13.33"},
        {"chr": "12", "start": 12600000, "end": 14600000, "band": "p13.1"},
        {"chr": "12", "start": 26300000, "end": 27600000, "band": "p11.23"},
        {"chr": "12", "start": 30500000, "end": 33200000, "band": "p11.21"},
        {"chr": "X", "start": 0, "end": 156040895, "band": "p11"},
        {"chr": "Y", "start": 0, "end": 57227415, "band": "p11"},
    ]


def run_tests():
    """Run built-in checks for parser, output header, coordinate sorting, and review routing."""
    raw_rows = [
        {"Sample_ID": "DEL", "ISCN": "46,XY,del(5)(q13q33)", "Ignored": "extra"},
        {"Sample_ID": "DUP", "ISCN": "46,XX,dup(7)(q21q31.33)"},
        {"Sample_ID": "GAIN", "ISCN": "47,XY,+8"},
        {"Sample_ID": "LOSS", "ISCN": "45,XX,-10"},
        {"Sample_ID": "MULTI", "ISCN": "46,XX,del(6)(q21q25.3),dup(6)(p25.3p22.1)"},
        {"Sample_ID": "SORTED", "ISCN": "46,XX,del(5)(q33q13)"},
        {"Sample_ID": "DECIMAL_PREFIX", "ISCN": "46,XX,dup(12)(p13p11.2)"},
        {"Sample_ID": "BALANCED", "ISCN": "46,XX,t(9;22)(q34;q11.2)"},
        {"Sample_ID": "AMBIG", "ISCN": "46,XX,del(5)(q13?q33)"},
        {"Sample_ID": "MIXED", "ISCN": "46,XX,t(9;22)(q34;q11.2),dup(7)(q21q31.33)"},
    ]
    plink, review, _debug = convert(raw_rows, _fake_cytobands())

    assert PLINK_COLUMNS == ["FID", "IID", "CHR", "BP1", "BP2", "TYPE (# of copies)", "SCORE", "SITES"]
    assert any(r["FID"] == "DEL" and r[TYPE_COLUMN] == 1 for r in plink)
    assert any(r["FID"] == "DUP" and r[TYPE_COLUMN] == 3 for r in plink)
    assert any(r["FID"] == "GAIN" and r[TYPE_COLUMN] == 3 for r in plink)
    assert any(r["FID"] == "LOSS" and r[TYPE_COLUMN] == 1 for r in plink)
    assert sum(1 for r in plink if r["FID"] == "MULTI") == 2
    assert any(r["FID"] == "MULTI" and r["CHR"] == "6" and r["BP1"] == 105000000 and r["BP2"] == 160600000 and r[TYPE_COLUMN] == 1 for r in plink)
    assert any(r["FID"] == "MULTI" and r["CHR"] == "6" and r["BP1"] == 0 and r["BP2"] == 30500000 and r[TYPE_COLUMN] == 3 for r in plink)
    assert any(r["FID"] == "DECIMAL_PREFIX" and r["BP1"] == 0 and r["BP2"] == 33200000 for r in plink)
    assert all(isinstance(r["BP1"], int) and isinstance(r["BP2"], int) and r["BP1"] < r["BP2"] for r in plink)
    assert all(r["FID"] == r["IID"] and r["FID"] for r in plink)
    validate_plink_rows(plink)
    assert not any(r["FID"] in {"BALANCED", "AMBIG"} for r in plink)
    assert any(r["Sample_ID"] == "BALANCED" and r["Status"] in {"FLAGGED", "SKIPPED"} for r in review)
    assert any(r["Sample_ID"] == "AMBIG" and r["Status"] == "FLAGGED" for r in review)
    assert any(r["FID"] == "MIXED" and r[TYPE_COLUMN] == 3 for r in plink)
    assert any(r["Sample_ID"] == "MIXED" and r["Reason"] == "unsupported_or_balanced_event" for r in review)

    with tempfile.TemporaryDirectory() as tmpdir:
        write_outputs(tmpdir, plink, review, _debug)
        backend = validate_written_plink_tsv(Path(tmpdir) / "plink_cnv_output.tsv")
        assert backend in {"pandas", "csv"}
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
        validate_written_plink_tsv(Path(a.outdir) / "plink_cnv_output.tsv")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    dels = sum(1 for r in plink if r[TYPE_COLUMN] == 1)
    dups = sum(1 for r in plink if r[TYPE_COLUMN] == 3)
    whole = sum(1 for r in debug if r["Event_Type"] in {"whole_gain", "whole_loss"} and r["Status"] == "EMITTED")
    print(f"Input samples: {len({r['Sample_ID'] for r in raw})}")
    print(f"CNV rows emitted: {len(plink)}")
    print(f"Flagged/skipped event records: {sum(1 for r in review if r['Status'] != 'EMITTED')}")
    print(f"Deletion rows: {dels}")
    print(f"Duplication rows: {dups}")
    print(f"Whole-chromosome gain/loss rows: {whole}")


if __name__ == "__main__":
    main()
