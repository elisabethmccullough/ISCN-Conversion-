# ISCN to PLINK CNV Converter

This repository contains a small command-line script that converts safe ISCN-like copy-number events into a PLINK-style CNV table.

You **do not need to manually merge** your raw ISCN document with the UCSC cytoband document. Keep them as two separate files and pass both paths to the script:

1. your raw sample file with `Sample_ID` and `ISCN` columns; and
2. the UCSC `cytoBand.txt.gz` reference file or URL.

The script reads both files, performs the cytoband lookup, and writes the merged/conversion outputs for you.


## Which file is the main raw file?

The **main raw file** is the file that contains your own samples and their ISCN strings. In this repository, `sample_raw.tsv` is only a tiny example file so you can test the script before using your real data.

When you run the script, the file after `--raw` is the main raw file. For example:

```powershell
py .\iscn_to_plink_cnv.py --raw .\sample_raw.tsv --cytoband .\sample_cytoband.tsv --outdir .\out
```

In that example, `sample_raw.tsv` is the raw file. For your real analysis, replace it with the path to your own CSV or TSV file, such as:

```powershell
py .\iscn_to_plink_cnv.py --raw .\my_patient_iscn_file.tsv --cytoband .\cytoBand.txt.gz --outdir .\results
```

Your raw file should look like this, with these exact column names in the first row:

```text
Sample_ID	ISCN
SIM001	46,XY,del(5)(q13q33)
SIM002	46,XX,dup(7)(q21q31.33)
```

Do **not** use the cytoband file as the raw file. The cytoband file goes after `--cytoband`; your sample/ISCN table goes after `--raw`.

## Required input files

### 1. Raw ISCN file

Create a CSV or TSV file with at least these two column names exactly:

```text
Sample_ID	ISCN
SIM001	46,XY,del(5)(q13q33)
SIM002	46,XX,dup(7)(q21q31.33)
```

Notes:

- In Excel, save as **CSV UTF-8** or **Text (Tab delimited)**.
- Keep the column names as `Sample_ID` and `ISCN`.
- Do not paste the cytoband table into this file.

### 2. UCSC cytoband file

Use the hg38 UCSC cytoband file:

```text
https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
```

If your network blocks direct URL access from Python, download the file first and use the local file path instead.

## Windows PowerShell quick start

Open **PowerShell** in the folder where you want to work, then run the commands below.

### Step 1: Check that Python is available

```powershell
py --version
```

If that fails, install Python from <https://www.python.org/downloads/windows/> and check the box that says **Add python.exe to PATH** during installation.

### Step 2: Go to the project folder

Replace the path below with the folder where you saved this repository.

```powershell
cd "C:\Users\YOUR_NAME\Documents\ISCN-Conversion-"
```

### Step 3: Optional, run the built-in smoke test

```powershell
py .\iscn_to_plink_cnv.py --run-tests
```

You should see:

```text
All tests passed.
```

### Step 4A: Run with the UCSC URL directly

Replace `my_raw_iscn.tsv` with your real raw file name.

```powershell
py .\iscn_to_plink_cnv.py `
  --raw .\my_raw_iscn.tsv `
  --cytoband "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz" `
  --outdir .\results
```

PowerShell uses the backtick character (`` ` ``) for line continuation. If you prefer a single-line command, use:

```powershell
py .\iscn_to_plink_cnv.py --raw .\my_raw_iscn.tsv --cytoband "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz" --outdir .\results
```

### Step 4B: If the URL is blocked, download cytoband first

```powershell
Invoke-WebRequest `
  -Uri "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz" `
  -OutFile ".\cytoBand.txt.gz"
```

Then run the converter with the local cytoband file:

```powershell
py .\iscn_to_plink_cnv.py --raw .\my_raw_iscn.tsv --cytoband .\cytoBand.txt.gz --outdir .\results
```

## Output files

The output directory will contain three files:

| File | Purpose |
| --- | --- |
| `plink_cnv_output.tsv` | Safe converted CNV rows only. Columns are exactly `FID`, `IID`, `CHR`, `BP1`, `BP2`, `TYPE (# of copies)`, `SCORE`, `SITES`. |
| `plink_cnv_output.csv` | Convenience comma-delimited copy of the main PLINK-style output for downloading/opening in spreadsheet software. |
| `review_flags.tsv` | Every input row/event that was emitted, approximately emitted, skipped, or flagged, with a reason. |
| `optional_debug_extracted_events.tsv` | Event-level debug table before final conversion, including `TYPE (# of copies)`, `Confidence`, and `Notes`. |


The main `plink_cnv_output.tsv` is a plain tab-delimited text file. You can save or rename it as `.txt` without changing the file contents. It should not contain Excel formulas, spreadsheet error values such as `#NAME?` or `#REF!`, Git patch text, markdown, or notes.

Some ambiguous-but-still-mappable ISCN events are emitted to the main PLINK table with a review/debug trail. For example, a deletion like `del(5)(q?q33)` can be emitted approximately from the chromosome 5 q-arm start through the chromosome end when the cytoband reference supports that inference, and a duplication like `dup(7)(q21q31~q32)` can be emitted through the broader q32 endpoint. These records are marked as `EMITTED_APPROXIMATE` in `review_flags.tsv` and `optional_debug_extracted_events.tsv`.

After a successful run, PowerShell will print summary counts like:

```text
Input samples: 7
CNV rows emitted: 8
Flagged/skipped event records: 2
Deletion rows: 4
Duplication rows: 4
Whole-chromosome gain/loss rows: 2
```

## Viewing the results in PowerShell

```powershell
Get-Content .\results\plink_cnv_output.tsv
Get-Content .\results\review_flags.tsv
```

You can also open the TSV files in Excel. If Excel does not split columns automatically, use **Data > From Text/CSV** and choose tab as the delimiter.

## Common errors and next steps

### `Missing required arguments: --raw, --cytoband`

You ran the script without telling it where your raw file and cytoband file are. Re-run with both required arguments:

```powershell
py .\iscn_to_plink_cnv.py --raw .\my_raw_iscn.tsv --cytoband .\cytoBand.txt.gz --outdir .\results
```

### `Raw file must include columns: Sample_ID and ISCN`

Your raw file does not have the required column headers. Open it in Excel or a text editor and make sure the first row contains exactly:

```text
Sample_ID	ISCN
```

For CSV files, the first row should be:

```text
Sample_ID,ISCN
```

### URL or 403/network errors

Your computer or institution may block Python from downloading the UCSC file directly. Download the cytoband file in a browser or with `Invoke-WebRequest`, then pass the local `.gz` file to `--cytoband`.

### No rows appear in `plink_cnv_output.tsv`

Check `review_flags.tsv`. The script is intentionally conservative and will flag ambiguous or unsupported events instead of guessing.

## Example using the included sample files

```powershell
py .\iscn_to_plink_cnv.py --raw .\sample_raw.tsv --cytoband .\sample_cytoband.tsv --outdir .\out
Get-Content .\out\plink_cnv_output.tsv
```
