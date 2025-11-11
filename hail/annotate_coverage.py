#!/usr/bin/env python3
import argparse
import logging
import math
import os
import re
import hail as hl
from os.path import dirname, abspath
from hail.utils.java import info

# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("annotate_coverage")
logger.setLevel(logging.INFO)


# -----------------------------------------------------------------------------
# Multi-way join helper
# -----------------------------------------------------------------------------
def multi_way_union_mts(mts: list, tmp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """Join multiple Hail MatrixTables progressively with checkpointing."""
    if not mts:
        raise ValueError("No MatrixTables provided to multi_way_union_mts()")

    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        n_jobs = math.ceil(len(staging) / chunk_size)
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []
        for i in range(n_jobs):
            to_merge = staging[chunk_size * i: chunk_size * (i + 1)]
            info(f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs")

            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda _: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            out_stage = os.path.join(tmp_dir, f"stage_{stage}_job_{i}.ht")
            info(f"Checkpointing to: {out_stage}")
            next_stage.append(merged.checkpoint(out_stage, overwrite=True))

        info(f"done stage {stage}")
        stage += 1
        staging = next_stage

    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
def main(args):
    # ---- 1. Prepare directories ----
    cwd = os.getcwd()
    tmp_dir = os.path.join(cwd, "temp")
    os.makedirs(tmp_dir, exist_ok=True)

    logger.info(f"[INIT] Working directory: {cwd}")
    logger.info(f"[INIT] Temp directory: {tmp_dir}")

    # ---- 2. Initialize Hail ----
    try:
        hl.init(
            tmp_dir=f"file:///{tmp_dir}",  # ✅ triple slashes
            local_tmpdir=tmp_dir,
            log=os.path.join(tmp_dir, "hail.log")
        )
        logger.info(f"[INIT] Hail initialized successfully. Version: {hl.__version__}")
    except Exception as e:
        logger.error(f"[FATAL] Hail init failed: {e}")
        raise

    # ---- 3. Load arguments ----
    input_tsv = abspath(args.input_tsv)
    output_ht = abspath(args.output_ht)
    chunk_size = args.chunk_size
    overwrite = args.overwrite

    logger.info(f"[ARGS] Input list: {input_tsv}")
    logger.info(f"[ARGS] Output HT: {output_ht}")
    logger.info(f"[ARGS] Chunk size: {chunk_size}")
    logger.info(f"[ARGS] Overwrite: {overwrite}")

    # ---- 4. Read TSV ----
    mt_list = []
    with open(input_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            items = line.split("\t")
            sample, cov_path = items[0:2]
            if not os.path.exists(cov_path):
                logger.warning(f"[SKIP] Coverage file missing: {cov_path}")
                continue
            logger.info(f"[LOAD] Importing coverage file: {cov_path}")

            mt = hl.import_matrix_table(
                cov_path,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
            ).drop("target")

            mt = mt.rename({"x": "coverage"})
            mt = mt.key_cols_by(s=sample)
            mt_list.append(mt)

    if not mt_list:
        raise RuntimeError("[ERROR] No valid MatrixTables loaded!")

    # ---- 5. Join coverage MTs ----
    out_dir = dirname(output_ht)
    temp_out_dir = os.path.join(out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)
    temp_out_dir_hail = f"file:///{temp_out_dir}"

    logger.info(f"[JOIN] Using Hail temp dir: {temp_out_dir_hail}")
    cov_mt = multi_way_union_mts(mt_list, temp_out_dir_hail, chunk_size)
    n_samples = cov_mt.count_cols()
    logger.info(f"[JOIN] Combined {n_samples} samples successfully")

    # ---- 6. Annotate coverage ----
    logger.info("[ANNOTATE] Adding row-level coverage annotations...")
    cov_mt = cov_mt.annotate_rows(
        locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float(hl.agg.count_where(cov_mt.coverage > 100) / n_samples),
        over_1000=hl.float(hl.agg.count_where(cov_mt.coverage > 1000) / n_samples),
    )

    # ---- 7. Write outputs ----
    output_mt = re.sub(r"\.ht$", ".mt", output_ht)
    output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
    output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

    logger.info("[WRITE] Exporting sample-level coverage...")
    sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
    sample_mt.coverage.export(output_samples)

    logger.info("[WRITE] Writing coverage matrix and table...")
    cov_mt.write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.rows()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    cov_ht.export(output_tsv)

    logger.info("[DONE] All files successfully written.")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine individual per-base coverage files into Hail MT/HT with annotations."
    )
    parser.add_argument("-i", "--input_tsv", required=True, help="Input TSV: sample_id<TAB>coverage_file")
    parser.add_argument("-o", "--output_ht", required=True, help="Output .ht path")
    parser.add_argument("--chunk_size", type=int, default=100, help="Join chunk size")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()
    main(args)
