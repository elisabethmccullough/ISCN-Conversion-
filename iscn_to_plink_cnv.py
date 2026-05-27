#!/usr/bin/env python3
"""Conservative ISCN-to-PLINK CNV converter.

Reads a raw file with Sample_ID and ISCN columns, plus a UCSC cytoband table,
extracts safe copy-number events, maps to base-pair coordinates, and writes:
- plink_cnv_output.tsv
- review_flags.tsv
- optional_debug_extracted_events.tsv
"""

from __future__ import annotations

import argparse
import io
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


PLINK_COLUMNS = ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"]


def _read_delimited(path_or_url: str) -> pd.DataFrame:
    """Read CSV/TSV/whitespace-delimited text from path or URL."""
    src = path_or_url.strip()
    if src.startswith(("http://", "https://")):
        content = pd.read_csv(src, sep=None, engine="python", compression="infer")
        return content
    try:
        return pd.read_csv(src, sep=None, engine="python")
    except Exception:
        return pd.read_csv(src, delim_whitespace=True)


def load_raw_data(path_or_url: str) -> pd.DataFrame:
    """Load raw ISCN input and validate required columns."""
    df = _read_delimited(path_or_url)
    required = {"Sample_ID", "ISCN"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required raw columns: {sorted(missing)}")
    return df[["Sample_ID", "ISCN"]].copy()


def load_cytoband(path_or_url: str) -> pd.DataFrame:
    """Load cytoband table and normalize columns."""
    try:
        df = pd.read_csv(path_or_url, sep="\t", header=None, compression="infer")
    except Exception:
        df = _read_delimited(path_or_url)

    if len(df.columns) < 4:
        raise ValueError("Cytoband file must contain at least 4 columns")

    df = df.iloc[:, :5].copy()
    cols = ["chrom", "start", "end", "band", "stain"][: len(df.columns)]
    df.columns = cols
    if "stain" not in df.columns:
        df["stain"] = ""

    df["chrom"] = df["chrom"].astype(str)
    df["band"] = df["band"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end", "chrom", "band"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["chr"] = df["chrom"].str.replace("^chr", "", regex=True)
    return df


def parse_del_dup_events(iscn: str) -> List[Dict[str, str]]:
    """Extract simple del/dup events of form del(chr)(bandStartbandEnd)."""
    events = []
    patt = re.compile(r"\b(del|dup)\(([^()]+)\)\(([pq][0-9.]+)([pq][0-9.]+)\)")
    for m in patt.finditer(iscn):
        etype, chrom, b1, b2 = m.groups()
        events.append(
            {
                "extracted_event": m.group(0),
                "event_type": etype,
                "chr": chrom.strip(),
                "band_start": b1.strip(),
                "band_end": b2.strip(),
            }
        )
    return events


def parse_whole_chromosome_events(iscn: str) -> List[Dict[str, str]]:
    """Extract whole-chromosome gains/losses like +8, -7, +X, -Y."""
    events = []
    for m in re.finditer(r"(?<![0-9A-Za-z])([+-])(\d+|X|Y)(?![0-9A-Za-z.])", iscn):
        sign, chrom = m.groups()
        events.append(
            {
                "extracted_event": m.group(0),
                "event_type": "whole_gain" if sign == "+" else "whole_loss",
                "chr": chrom,
                "band_start": "",
                "band_end": "",
            }
        )
    return events


def lookup_cytoband_coordinates(cyto: pd.DataFrame, chr_code: str, band: str) -> Tuple[Optional[int], Optional[int], str]:
    """Lookup band coordinates with exact first, then safe prefix expansion."""
    chr_code = str(chr_code).replace("chr", "")
    cdf = cyto[cyto["chr"] == chr_code]
    if cdf.empty:
        return None, None, "chromosome_not_found"

    exact = cdf[cdf["band"] == band]
    if not exact.empty:
        return int(exact["start"].min()), int(exact["end"].max()), "high"

    pref = cdf[cdf["band"].str.startswith(f"{band}.", na=False)]
    if pref.empty:
        return None, None, "band_not_found"
    return int(pref["start"].min()), int(pref["end"].max()), "medium"


def validate_output_row(row: Dict) -> Tuple[bool, str]:
    """Validate PLINK output row safety constraints."""
    try:
        bp1 = int(row["BP1"])
        bp2 = int(row["BP2"])
    except Exception:
        return False, "non_integer_coordinates"
    if bp1 >= bp2:
        return False, "reversed_or_zero_length_coordinates"
    if str(row["CHR"]) in {"", "nan", "None"}:
        return False, "missing_chromosome"
    if int(row["TYPE"]) not in {0, 1, 3}:
        return False, "unsupported_type"
    return True, "ok"


def write_outputs(outdir: Path, plink_rows: List[Dict], review_rows: List[Dict], debug_rows: List[Dict]) -> None:
    """Write required and optional output tables."""
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(plink_rows, columns=PLINK_COLUMNS).to_csv(outdir / "plink_cnv_output.tsv", sep="\t", index=False)
    pd.DataFrame(review_rows).to_csv(outdir / "review_flags.tsv", sep="\t", index=False)
    pd.DataFrame(debug_rows).to_csv(outdir / "optional_debug_extracted_events.tsv", sep="\t", index=False)


def _contains_ambiguity(text: str) -> bool:
    return any(sym in text for sym in ["?", "~"]) or "q?q" in text.lower()


def convert(raw_df: pd.DataFrame, cyto_df: pd.DataFrame):
    plink_rows, review_rows, debug_rows = [], [], []
    for _, r in raw_df.iterrows():
        sid = str(r["Sample_ID"])
        raw = str(r["ISCN"])

        if _contains_ambiguity(raw):
            review_rows.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": "", "Status": "FLAGGED", "Reason": "ambiguous_symbols", "Confidence": "low", "Notes": "Contains ? or ~"})
            continue

        events = parse_del_dup_events(raw) + parse_whole_chromosome_events(raw)

        if not events:
            reason = "non_cnv_or_no_safe_event"
            if re.search(r"\b(t\(|add\(|mar\b|der\()", raw):
                reason = "unsupported_or_balanced_event"
            review_rows.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": "", "Status": "SKIPPED", "Reason": reason, "Confidence": "low", "Notes": "No safely convertible CNV event"})
            continue

        for ev in events:
            chr_code = ev["chr"]
            event_type = ev["event_type"]
            extracted = ev["extracted_event"]
            status, reason, conf, notes = "SKIPPED", "", "low", ""
            bp1 = bp2 = cn_type = None

            if event_type in {"del", "dup"}:
                s1, e1, c1 = lookup_cytoband_coordinates(cyto_df, chr_code, ev["band_start"])
                s2, e2, c2 = lookup_cytoband_coordinates(cyto_df, chr_code, ev["band_end"])
                if None in {s1, e1, s2, e2}:
                    reason = "cytoband_lookup_failed"
                else:
                    bp1, bp2 = int(min(s1, s2)), int(max(e1, e2))
                    cn_type = 1 if event_type == "del" else 3
                    conf = "medium" if "medium" in {c1, c2} else "high"
            else:
                cdf = cyto_df[cyto_df["chr"] == str(chr_code)]
                if cdf.empty:
                    reason = "chromosome_not_found"
                else:
                    bp1, bp2 = int(cdf["start"].min()), int(cdf["end"].max())
                    cn_type = 3 if event_type == "whole_gain" else 1
                    conf = "high"

            if bp1 is not None and bp2 is not None and cn_type is not None:
                row = {"FID": sid, "IID": sid, "CHR": str(chr_code), "BP1": bp1, "BP2": bp2, "TYPE": cn_type, "SCORE": 999, "SITES": 999}
                ok, vreason = validate_output_row(row)
                if ok:
                    plink_rows.append(row)
                    status, reason = "EMITTED", "safe_conversion"
                else:
                    status, reason = "FLAGGED", vreason
            else:
                status = "FLAGGED"
                reason = reason or "unsupported_event"

            review_rows.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": extracted, "Status": status, "Reason": reason, "Confidence": conf, "Notes": notes})
            debug_rows.append({"Sample_ID": sid, "Raw_ISCN": raw, "Extracted_Event": extracted, "Event_Type": event_type, "CHR": chr_code, "Band_Start": ev["band_start"], "Band_End": ev["band_end"], "BP1": bp1, "BP2": bp2, "TYPE": cn_type, "Status": status, "Reason": reason})

    return plink_rows, review_rows, debug_rows


def run_tests() -> None:
    """Run built-in sanity tests with synthetic data."""
    cyto = pd.DataFrame(
        [
            ["chr5", 67000000, 70000000, "q13", "gneg"],
            ["chr5", 70000000, 100000000, "q14", "gpos"],
            ["chr5", 150000000, 160500000, "q33", "gneg"],
            ["chr6", 0, 2300000, "p25.3", "gneg"],
            ["chr6", 25000000, 30500000, "p22.1", "gpos"],
            ["chr6", 105000000, 110000000, "q21", "gneg"],
            ["chr6", 158000000, 160600000, "q25.3", "gneg"],
            ["chr7", 77900000, 80000000, "q21", "gneg"],
            ["chr7", 120000000, 127500000, "q31.33", "gneg"],
            ["chr8", 0, 145000000, "p11", "gneg"],
            ["chr10", 0, 133000000, "p11", "gneg"],
            ["chrY", 0, 57200000, "p11", "gneg"],
        ],
        columns=["chrom", "start", "end", "band", "stain"],
    )
    cyto["chr"] = cyto["chrom"].str.replace("^chr", "", regex=True)

    raw = pd.DataFrame(
        {
            "Sample_ID": ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"],
            "ISCN": [
                "46,XY,del(5)(q13q33)",
                "46,XX,dup(7)(q21q31.33)",
                "46,XX,del(6)(q21q25.3),dup(6)(p25.3p22.1)",
                "47,XY,+8",
                "45,XX,-10",
                "46,XX",
                "46,XX,t(9;22)(q34;q11.2)",
                "46,XX,del(5)(q13?q33)",
                "46,XX,del(99)(q13q33)",
                "46,XX,del(5)(q33q13)",
            ],
        }
    )

    plink, review, _ = convert(raw, cyto)
    out = pd.DataFrame(plink)

    assert ((out["FID"] == "S1") & (out["TYPE"] == 1)).any()
    assert ((out["FID"] == "S2") & (out["TYPE"] == 3)).any()
    assert (out["FID"] == "S3").sum() == 2
    assert ((out["FID"] == "S4") & (out["TYPE"] == 3)).any()
    assert ((out["FID"] == "S5") & (out["TYPE"] == 1)).any()
    assert not (out["FID"] == "S6").any()
    assert any(r["Sample_ID"] == "S7" and r["Status"] in {"SKIPPED", "FLAGGED"} for r in review)
    assert any(r["Sample_ID"] == "S8" and r["Status"] == "FLAGGED" for r in review)
    assert any(r["Sample_ID"] == "S9" and r["Status"] == "FLAGGED" for r in review)
    assert any(r["Sample_ID"] == "S10" and r["Status"] == "FLAGGED" for r in review)
    print("All tests passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Conservative ISCN to PLINK CNV converter")
    parser.add_argument("--raw", help="Path or URL to raw ISCN file")
    parser.add_argument("--cytoband", help="Path or URL to cytoband file")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--run-tests", action="store_true", help="Run built-in tests and exit")
    args = parser.parse_args()

    if args.run_tests:
        run_tests()
        return

    if not args.raw or not args.cytoband:
        raise SystemExit("--raw and --cytoband are required unless --run-tests is used")

    raw_df = load_raw_data(args.raw)
    cyto_df = load_cytoband(args.cytoband)
    plink_rows, review_rows, debug_rows = convert(raw_df, cyto_df)
    write_outputs(Path(args.outdir), plink_rows, review_rows, debug_rows)

    plink_df = pd.DataFrame(plink_rows)
    del_count = int((plink_df["TYPE"] == 1).sum()) if not plink_df.empty else 0
    dup_count = int((plink_df["TYPE"] == 3).sum()) if not plink_df.empty else 0
    whole_count = int(sum(1 for r in debug_rows if r["Event_Type"] in {"whole_gain", "whole_loss"} and r["Status"] == "EMITTED"))

    print(f"Input samples: {raw_df['Sample_ID'].nunique()}")
    print(f"CNV rows emitted: {len(plink_rows)}")
    print(f"Flagged/skipped event records: {sum(1 for r in review_rows if r['Status'] != 'EMITTED')}")
    print(f"Deletion rows: {del_count}")
    print(f"Duplication rows: {dup_count}")
    print(f"Whole-chromosome gain/loss rows: {whole_count}")


if __name__ == "__main__":
    main()
