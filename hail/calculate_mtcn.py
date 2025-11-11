#!/usr/bin/env python3
"""
calculate_mtcn.py
=================
Pipeline steps:
  1) Run Picard CollectWgsMetrics on a single CRAM to obtain nuclear coverage.
  2) Parse nuclear coverage (autosomes 1–22 and sex chromosomes X/Y).
  3) Compute mitochondrial coverage (mean and median) from a per-base TSV.
  4) Compute mtDNA copy number (mtCN):
        - mtCN_mean_mean     = 2 × mean_mt / mean_nuclear
        - mtCN_median_median = 2 × median_mt / median_nuclear
        - mtCN_final         = 2 × mean_mt / median_nuclear  (used for QC)
  5) Apply QC filters (based on mtCN_final) and write mtCN_summary.txt.
"""

import os
import io
import re
import argparse
import subprocess
import pandas as pd
import numpy as np


# -------------------------------------------------------------------------
# Step 1: Generate WGS metrics with Picard
# -------------------------------------------------------------------------
def generate_wgs_metrics(
    cram: str,
    ref_fasta: str,
    output_dir: str,
    picard_path: str,
    intervals: str | None = None,
    use_fast_algorithm: bool = True,
    read_length: int = 151,
) -> str:
    """
    Run Picard CollectWgsMetrics for a single CRAM.

    Outputs
    -------
    <output_dir>/<sample>.wgsMetrics.txt

    Notes
    -----
    - If `intervals` is provided, Picard computes on that subset only (e.g., autosomes + XY).
    - If `use_fast_algorithm` is True, passes USE_FAST_ALGORITHM=true to Picard.
    """
    os.makedirs(output_dir, exist_ok=True)
    sample = os.path.basename(cram).removesuffix(".cram")
    out_metrics = os.path.join(output_dir, f"{sample}.wgsMetrics.txt")

    cmd = [
        "java", "-Xmx8g", "-jar", picard_path, "CollectWgsMetrics",
        "VALIDATION_STRINGENCY=SILENT",
        f"INPUT={cram}",
        f"OUTPUT={out_metrics}",
        f"REFERENCE_SEQUENCE={ref_fasta}",
        "MINIMUM_MAPPING_QUALITY=20",
        "MINIMUM_BASE_QUALITY=20",
        f"READ_LENGTH={int(read_length)}",
    ]
    if intervals:
        cmd.append(f"INTERVALS={intervals}")
    if use_fast_algorithm:
        cmd.append("USE_FAST_ALGORITHM=true")

    print(f"[CMD] {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"[OK] WGS metrics generated for {sample}")
    except subprocess.CalledProcessError as e:
        print(f"[ERR] Picard failed for {sample}: {e}")
        raise RuntimeError(f"Picard CollectWgsMetrics failed for {sample}") from e

    return out_metrics


# -------------------------------------------------------------------------
# Step 2: Parse nuclear coverage (autosomes + XY)
# -------------------------------------------------------------------------
def parse_nuclear_coverage(metrics_file: str) -> tuple[str, float, float]:
    """
    Parse the first WgsMetrics table block from Picard output (exclude histograms).

    Returns
    -------
    (sample_name, nuclear_median_coverage, nuclear_mean_coverage)
    """
    with open(metrics_file, "r") as fh:
        text = fh.read()

    # Locate the first WgsMetrics table header and take the first data row
    m = re.search(r"## METRICS CLASS[^\n]*WgsMetrics[^\n]*\n", text, flags=re.MULTILINE)
    if not m:
        raise RuntimeError(f"Cannot locate WgsMetrics table in: {metrics_file}")

    block = text[m.end():]
    block = block.split("## HISTOGRAM", 1)[0].strip()
    lines = [ln for ln in block.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"WgsMetrics table has no data row in: {metrics_file}")

    header, first = lines[0], lines[1]
    df = pd.read_csv(io.StringIO(header + "\n" + first + "\n"), sep="\t", engine="python")

    def pick(colnames, candidates):
        for c in candidates:
            if c in colnames:
                return c
        raise KeyError(f"None of columns {candidates} found in {list(colnames)}")

    median_col = pick(df.columns, ["MEDIAN_COVERAGE"])
    mean_col   = pick(df.columns, ["MEAN_COVERAGE"])

    nuc_median = float(pd.to_numeric(df[median_col], errors="coerce").iloc[0])
    nuc_mean   = float(pd.to_numeric(df[mean_col],   errors="coerce").iloc[0])

    sample = os.path.basename(metrics_file).split(".")[0]
    return sample, nuc_median, nuc_mean


# -------------------------------------------------------------------------
# Step 3: Compute mitochondrial coverage
# -------------------------------------------------------------------------
def calculate_mt_coverage(mt_coverage_file: str) -> tuple[str, float, float]:
    """
    Compute mean and median mtDNA coverage from a per-base coverage TSV.

    Expected columns
    ----------------
    'coverage'
    """
    if not os.path.exists(mt_coverage_file):
        raise FileNotFoundError(f"MT coverage file not found: {mt_coverage_file}")

    df = pd.read_csv(mt_coverage_file, sep="\t")
    if "coverage" not in df.columns:
        raise ValueError("Expected 'coverage' column in mt coverage file.")

    cov = pd.to_numeric(df["coverage"], errors="coerce")
    cov_mean = float(cov.mean())
    cov_median = float(cov.median())

    sample = os.path.basename(mt_coverage_file).split(".")[0]
    print(f"[OK] {sample}: mtDNA mean={cov_mean:.2f}, median={cov_median:.2f}")
    return sample, cov_mean, cov_median


# -------------------------------------------------------------------------
# Step 4: Compute mtCN and apply QC filters
# -------------------------------------------------------------------------
def compute_mtcn(
    sample: str,
    mt_mean: float,
    mt_median: float,
    nuc_mean: float,
    nuc_median: float,
    args: argparse.Namespace,
) -> dict:
    """
    Compute mtDNA copy number and determine pass/fail.

    Formulas
    --------
    mtCN_mean_mean     = 2 × mean_mt / mean_nuclear
    mtCN_median_median = 2 × median_mt / median_nuclear
    mtCN_final         = 2 × mean_mt / median_nuclear  (used for QC)
    """
    mtcn_mean_mean = 2 * mt_mean / nuc_mean if (nuc_mean and np.isfinite(nuc_mean) and nuc_mean > 0) else np.nan
    mtcn_median_median = 2 * mt_median / nuc_median if (nuc_median and np.isfinite(nuc_median) and nuc_median > 0) else np.nan
    mtcn_final = 2 * mt_mean / nuc_median if (nuc_median and np.isfinite(nuc_median) and nuc_median > 0) else np.nan

    # Placeholder: contamination may be added by upstream logic if available
    contam = 0.0

    pass_filter = (
        (mt_mean >= args.min_mt_cov)
        and (np.isfinite(mtcn_final))
        and (args.min_mtcn <= mtcn_final <= args.max_mtcn)
        and (contam <= args.max_contam)
    )

    print(f"[RESULT] {sample}: mtCN_final={mtcn_final:.2f}, pass={pass_filter}")
    return {
        "Sample_ID": sample,
        "Mean_mtDNA_Coverage": mt_mean,
        "Median_mtDNA_Coverage": mt_median,
        "Mean_Nuclear_Coverage": nuc_mean,
        "Median_Nuclear_Coverage": nuc_median,
        "mtCN_mean_mean": mtcn_mean_mean,
        "mtCN_median_median": mtcn_median_median,
        "mtCN_final": mtcn_final,
        "Contamination": contam,
        "Pass_Filter": pass_filter,
    }


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main(args: argparse.Namespace) -> None:
    os.makedirs(args.output, exist_ok=True)

    print("[INFO] Step 1 – Generating WGS metrics")
    metrics_file = generate_wgs_metrics(
        cram=args.cram,
        ref_fasta=args.ref_fasta,
        output_dir=args.output,
        picard_path=args.picard,
        intervals=args.intervals,
        use_fast_algorithm=args.use_fast_algorithm,
        read_length=args.read_length,
    )

    print("[INFO] Step 2 – Parsing nuclear coverage (autosomes + XY)")
    sample, nuc_median, nuc_mean = parse_nuclear_coverage(metrics_file)

    print("[INFO] Step 3 – Calculating mtDNA coverage")
    _, mt_mean, mt_median = calculate_mt_coverage(args.mt_coverage)

    print("[INFO] Step 4 – Computing mtCN and applying QC filters")
    record = compute_mtcn(sample, mt_mean, mt_median, nuc_mean, nuc_median, args)

    out_tsv = os.path.join(args.output, "mtCN_summary.txt")
    pd.DataFrame([record]).to_csv(out_tsv, sep="\t", index=False)
    print(f"[DONE] Summary written to {out_tsv}")
    print(f"[SUMMARY] PASS = {record['Pass_Filter']}")


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate mtDNA copy number (mtCN) from a CRAM and per-base mt coverage TSV.\n"
            "Final mtCN = 2 × mean(mtDNA coverage) / median(nuclear coverage).\n"
            "Tip: use --intervals autosome_XY.interval_list to speed up Picard and align to autosomes + X/Y."
        )
    )
    parser.add_argument("--cram", required=True, help="Input CRAM file.")
    parser.add_argument("--crai", required=True, help="Index CRAI file (not used here; kept for interface parity).")
    parser.add_argument("--ref_fasta", required=True, help="Reference genome FASTA.")
    parser.add_argument("--mt_coverage", required=True, help="Per-base mt coverage TSV with a 'coverage' column.")
    parser.add_argument("--output", required=True, help="Output directory.")
    parser.add_argument("--picard", required=True, help="Path to Picard JAR.")

    # Performance and scope
    parser.add_argument(
        "--intervals",
        default=None,
        help="Picard interval_list to restrict WGS metrics (e.g., autosomes + XY)."
    )
    parser.add_argument(
        "--use_fast_algorithm",
        action="store_true",
        help="Enable Picard USE_FAST_ALGORITHM=true (faster on high coverage)."
    )
    parser.add_argument("--read_length", type=int, default=151, help="Read length hint for Picard.")

    # QC thresholds
    parser.add_argument("--min_mt_cov", type=int, default=100, help="Minimum mtDNA mean coverage.")
    parser.add_argument("--min_mtcn", type=float, default=50, help="Minimum mtCN threshold.")
    parser.add_argument("--max_mtcn", type=float, default=500, help="Maximum mtCN threshold.")
    parser.add_argument("--max_contam", type=float, default=0.02, help="Maximum contamination allowed.")

    main(parser.parse_args())
