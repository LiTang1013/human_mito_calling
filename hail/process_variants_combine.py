#!/usr/bin/env python3

import os
import io
import json
import argparse
import pandas as pd


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #
def load_complex_db(path: str, key_cols: list[str], val_cols: list[str], delimiter: str = "\t") -> dict:
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#")
    missing = [c for c in key_cols + val_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: {os.path.basename(path)} missing columns: {missing}")
    idx = pd.MultiIndex.from_frame(df[key_cols])
    df["__values__"] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df["__values__"].values, index=idx).to_dict()


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description="Process a VEP-annotated multi-sample VCF and generate summary tables.")
    ap.add_argument("--vep-vcf", required=True, help="Path to one VEP-annotated VCF.")
    ap.add_argument("--final-output-dir", required=True, help="Directory for output tables.")

    # External databases (required)
    ap.add_argument("--gnomadcache", required=True)
    ap.add_argument("--clinvarcache", required=True)
    ap.add_argument("--mitomap-polycache", required=True)
    ap.add_argument("--mitomap-diseasecache", required=True)
    ap.add_argument("--helixcache", required=True)
    ap.add_argument("--mitimpactcache", required=True)
    ap.add_argument("--mitotipcache", required=True)
    ap.add_argument("--hmtvarcache", required=True)

    # Mode
    ap.add_argument("--pipeline-mode", choices=["population", "disease"], required=True)
    ap.add_argument(
        "--disease-meta-file",
        help="TSV with columns SampleID and Category; rows with Category=Proband are kept when --pipeline-mode=disease."
    )

    args = ap.parse_args()
    outdir = args.final_output_dir
    os.makedirs(outdir, exist_ok=True)

    # -------------------- disease-mode: select probands -------------------- #
    proband_set: set[str] | None = None
    if args.pipeline_mode == "disease":
        if not args.disease_meta_file:
            raise SystemExit("ERROR: --pipeline-mode=disease requires --disease-meta-file.")
        meta = pd.read_csv(args.disease_meta_file, sep="\t", dtype=str)
        if meta.empty:
            raise SystemExit("ERROR: disease_meta_file is empty.")
        meta.columns = [c.strip() for c in meta.columns]
        needed = {"SampleID", "Category"}
        if not needed.issubset(meta.columns):
            raise SystemExit("ERROR: disease_meta_file must contain columns: SampleID, Category.")
        proband_set = set(
            meta.loc[meta["Category"].astype(str).str.lower() == "proband", "SampleID"].astype(str)
        )
        if not proband_set:
            raise SystemExit("ERROR: No rows with Category=Proband in disease_meta_file.")

    # -------------------- external annotation DBs ------------------------- #
    gnomad_db          = load_complex_db(args.gnomadcache,          ["ref", "position", "alt"], ["max_observed_heteroplasmy", "AF_hom", "AF_het"])
    mitomap_poly_db    = load_complex_db(args.mitomap_polycache,    ["ref", "pos", "alt"],      ["gbcnt"])
    mitomap_disease_db = load_complex_db(args.mitomap_diseasecache, ["ref", "pos", "alt"],      ["status", "homoplasmy", "heteroplasmy", "disease"])
    apogee_db          = load_complex_db(args.mitimpactcache,       ["Ref", "Start", "Alt"],    ["APOGEE1", "APOGEE2"])
    hmtvar_db          = load_complex_db(args.hmtvarcache,          ["REF", "POS", "ALT"],      ["HmtVar"])

    # MitoTIP quartile → verbal class
    mitotip_map = {
        "Q1": "likely pathogenic",
        "Q2": "possibly pathogenic",
        "Q3": "possibly benign",
        "Q4": "likely benign",
    }
    mitotip_df = pd.read_csv(args.mitotipcache, sep="\t", dtype=str)
    missing = [c for c in ["Quartile", "rCRS", "Position", "Alt"] if c not in mitotip_df.columns]
    if missing:
        raise SystemExit(f"ERROR: mitotipcache missing columns: {missing}")
    mitotip_df["prediction"] = mitotip_df["Quartile"].map(mitotip_map)
    mitotip_db = pd.Series(
        mitotip_df["prediction"].values,
        index=list(zip(mitotip_df["rCRS"], mitotip_df["Position"], mitotip_df["Alt"]))
    ).to_dict()

    # HelixMTdb parsing
    helix_df = pd.read_csv(args.helixcache, sep="\t", dtype=str)
    need_helix = {"locus", "alleles", "AF_hom", "AF_het", "max_ARF", "ref", "alt", "pos"}
    # derive ref/alt/pos when necessary
    helix_df = helix_df[helix_df["alleles"].str.count(",") == 1].copy()
    helix_df["pos"] = helix_df["locus"].str.split("chrM:").str[1]
    a = helix_df["alleles"].str.split('"').str
    helix_df["ref"] = a[1]
    helix_df["alt"] = a[3]
    for col in ["AF_hom", "max_ARF"]:
        helix_df[col] = pd.to_numeric(helix_df[col], errors="coerce").fillna(0.0)
    helix_df["max_het"] = helix_df.apply(lambda r: 1.0 if r["AF_hom"] > 0 else r["max_ARF"], axis=1)
    helix_db = pd.Series(
        list(zip(helix_df["max_het"], helix_df["AF_hom"], helix_df.get("AF_het", 0.0))),
        index=pd.MultiIndex.from_frame(helix_df[["ref", "pos", "alt"]])
    ).to_dict()

    # ClinVar: SNVs only + exclude "Conflicting interpretations of pathogenicity"
    clinvar_df = pd.read_csv(args.clinvarcache, sep="\t", dtype=str)
    needed = {"GRCh38Location", "Canonical SPDI", "Germline classification"}
    if not needed.issubset(clinvar_df.columns):
        raise SystemExit(f"ERROR: clinvarcache missing columns: {sorted(needed - set(clinvar_df.columns))}")
    clinvar_df["pos"] = clinvar_df["GRCh38Location"]
    spdi = clinvar_df["Canonical SPDI"].str.split(":", expand=True)
    clinvar_df["ref"] = spdi[2]
    clinvar_df["alt"] = spdi[3]
    clinvar_df = clinvar_df[
        (clinvar_df["ref"].str.len() == 1) &
        (clinvar_df["alt"].str.len() == 1) &
        (clinvar_df["ref"] != clinvar_df["alt"]) &
        (clinvar_df["Germline classification"] != "Conflicting interpretations of pathogenicity")
    ].copy()
    clinvar_db = pd.Series(
        clinvar_df["Germline classification"].values,
        index=pd.MultiIndex.from_frame(clinvar_df[["ref", "pos", "alt"]])
    ).to_dict()

    # -------------------- load VEP VCF -------------------- #
    with open(args.vep_vcf, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    if not lines:
        raise SystemExit("ERROR: VCF has no header/data after removing '##' meta lines.")
    vcf_df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str).rename(columns={"#CHROM": "CHROM"})
    if vcf_df.empty:
        raise SystemExit("ERROR: Empty VCF after header removal.")
    if "INFO" not in vcf_df.columns:
        raise SystemExit("ERROR: VCF does not contain an INFO column.")

    # Identify sample columns (from column index >= 9)
    sample_cols = list(vcf_df.columns[9:])
    if args.pipeline_mode == "disease":
        # Keep sample columns whose header begins with any proband ID (tolerate suffixes)
        kept = []
        for col in sample_cols:
            sid_prefix = col.split(":")[0] if ":" in col else col
            if any(sid_prefix.startswith(p) for p in proband_set):  # type: ignore[arg-type]
                kept.append(col)
        sample_cols = kept
        if not sample_cols:
            raise SystemExit("ERROR: No proband columns found in the VCF for disease mode.")

    # Parse VEP consequence fields from INFO ("|"-delimited)
    vep_cols = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
        "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
        "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "HGVS_OFFSET",
    ]
    info_split = vcf_df["INFO"].str.split("|", expand=True)
    # Safe assignment when the number of columns may be shorter
    for i, col in enumerate(vep_cols, start=1):
        if i < info_split.shape[1]:
            vcf_df[col] = info_split.iloc[:, i]
        else:
            vcf_df[col] = pd.NA

    # Long format: one row per sample/variant
    stacks = []
    for sc in sample_cols:
        sub = vcf_df.copy()
        sub["SampleID"] = sc
        sub["SAMPLE_DATA"] = sub[sc]
        stacks.append(sub)
    long_df = pd.concat(stacks, ignore_index=True)

    # Remove common hard filters
    long_df = long_df[~long_df["FILTER"].str.contains(
        "base_qual|strand_bias|weak_evidence|blacklisted_site|contamination|position", na=False
    )].copy()

    # FORMAT parsing helpers
    fmt = long_df["FORMAT"].str.split(":", expand=True)

    def pick(tag: str):
        """Return the column extracted from SAMPLE_DATA for a given FORMAT tag."""
        try:
            idx = list(fmt.iloc[0].values).index(tag)
            return long_df["SAMPLE_DATA"].str.split(":", expand=True).iloc[:, idx]
        except Exception:
            return pd.NA

    long_df["GT"] = pick("GT")
    long_df["AD"] = pick("AD")
    long_df["AF"] = pd.to_numeric(pick("AF"), errors="coerce")
    long_df["DP"] = pick("DP")

    # Annotation lookups keyed by (REF, POS, ALT)
    var_keys = list(zip(long_df["REF"], long_df["POS"], long_df["ALT"]))

    # gnomAD
    gvals = [gnomad_db.get(k, ("0", "0", "0")) for k in var_keys]
    long_df["gnomad_max_hl"] = [x[0] for x in gvals]
    long_df["gnomad_af_hom"] = pd.to_numeric([x[1] for x in gvals], errors="coerce")
    long_df["gnomad_af_het"] = pd.to_numeric([x[2] for x in gvals], errors="coerce")

    # Helix
    hvals = [helix_db.get(k, (0.0, 0.0, 0.0)) for k in var_keys]
    long_df["helix_max_hl"] = [x[0] for x in hvals]
    long_df["helix_af_hom"] = pd.to_numeric([x[1] for x in hvals], errors="coerce")
    long_df["helix_af_het"] = pd.to_numeric([x[2] for x in hvals], errors="coerce")

    # ClinVar
    long_df["clinvar_interp"] = [clinvar_db.get(k, "") for k in var_keys]

    # MITOMAP (polymorphisms + disease)
    long_df["mitomap_gbcnt"] = [mitomap_poly_db.get(k, ("0",))[0] for k in var_keys]
    long_df["mitomap_af"] = pd.to_numeric(long_df["mitomap_gbcnt"], errors="coerce") / 61134.0

    md = [mitomap_disease_db.get(k, ("", "", "", "")) for k in var_keys]
    long_df["mitomap_status"] = [d[0] for d in md]
    long_df["mitomap_plasmy"] = [f"{d[1]}/{d[2]}" if d[1] or d[2] else "" for d in md]
    long_df["mitomap_disease"] = [d[3] for d in md]

    # APOGEE (MITImpact)
    long_df["apogee_class"] = [
        str(apogee_db.get(k, "")).strip("()").replace("'", "").replace(", ", "/") for k in var_keys
    ]

    # MitoTIP
    long_df["mitotip_class"] = [mitotip_db.get(k, "") for k in var_keys]

    # HmtVar (JSON in a single field)
    def get_hmtvar(k):
        val = hmtvar_db.get(k)
        if val and val[0]:
            try:
                return json.loads(val[0]).get("pathogenicity", "")
            except json.JSONDecodeError:
                return ""
        return ""

    long_df["hmtvar_class"] = [get_hmtvar(k) for k in var_keys]

    # In-cohort AC and per-sample counts
    long_df["variant_key"] = long_df["POS"].astype(str) + ":" + long_df["REF"] + ":" + long_df["ALT"]
    ac = long_df["variant_key"].value_counts().rename("in_cohort_AC").reset_index().rename(columns={"index": "variant_key"})
    long_df = long_df.merge(ac, on="variant_key", how="left")
    long_df["in_cohort_AC"] = long_df["in_cohort_AC"].fillna(0).astype(int)

    freq = long_df["SampleID"].value_counts().rename("Freq").reset_index().rename(columns={"index": "SampleID"})
    long_df = long_df.merge(freq, on="SampleID", how="left")

    # Columns to emit (prefilter)
    base_cols = [
        "SampleID", "variant_key", "CHROM", "POS", "REF", "ALT", "FILTER",
        "GT", "AD", "AF", "DP",
        "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "HGVSp", "Codons", "VARIANT_CLASS",
        "gnomad_max_hl", "gnomad_af_hom", "gnomad_af_het",
        "apogee_class", "mitotip_class", "hmtvar_class",
        "helix_max_hl", "helix_af_hom", "helix_af_het",
        "mitomap_gbcnt", "mitomap_af", "mitomap_status", "mitomap_plasmy", "mitomap_disease",
        "clinvar_interp",
        "in_cohort_AC", "Freq",
    ]
    for c in base_cols:
        if c not in long_df.columns:
            long_df[c] = pd.NA
    pre = long_df[base_cols].copy()

    # Filenames
    if args.pipeline_mode == "disease":
        pre_path = os.path.join(outdir, "Proband_variant_list_prefiltering.txt")
        final_path = os.path.join(outdir, "Proband_variant_list.txt")
    else:
        pre_path = os.path.join(outdir, "variant_list_prefiltering.txt")
        final_path = os.path.join(outdir, "variant_list.txt")

    # Prefilter output (rename AF to Heteroplasmy)
    pre_out = pre.rename(columns={"AF": "Heteroplasmy"}).copy()
    pre_out.to_csv(pre_path, sep="\t", index=False, na_rep="")
    print(f"[+] Pre-filtering table written: {pre_path}")

    # Numeric conversions for filtering
    for col in ["gnomad_af_hom", "helix_af_hom", "mitomap_af", "Heteroplasmy", "Freq", "in_cohort_AC"]:
        if col in pre_out.columns:
            pre_out[col] = pd.to_numeric(pre_out[col], errors="coerce")

    # Final filter (kept as in your logic)
    filt = pre_out[
        (pre_out["gnomad_af_hom"] < 0.01) &
        (pre_out["helix_af_hom"] < 0.01) &
        (pre_out["mitomap_af"] < 0.01) &
        (pre_out["Consequence"] != "synonymous_variant") &
        (pre_out["Freq"] < 100) &
        (pre_out["Heteroplasmy"] > 0.05) &
        (~pre_out["clinvar_interp"].isin(["Benign", "Likely benign"]))
    ].copy()

    # Optional sort: disease-present first, then rarer in population resources
    if not filt.empty:
        filt["__mitomap_has_disease__"] = filt["mitomap_disease"].apply(lambda x: 0 if pd.isna(x) or x == "" else 1)
        filt = (
            filt.sort_values(
                by=["__mitomap_has_disease__", "gnomad_af_hom", "helix_af_hom", "mitomap_af"],
                ascending=[False, True, True, True],
            )
            .drop(columns="__mitomap_has_disease__")
        )

    filt.to_csv(final_path, sep="\t", index=False, na_rep="")
    print(f"[+] Final filtered table written: {final_path}")


if __name__ == "__main__":
    main()
