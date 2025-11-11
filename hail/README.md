# Hail Script Toolkit
This directory contains lightweight Python utilities for post-calling mitochondrial variant analysis. The tools are designed to run after variant calling step and cover the typical downstream flow:

### mtCN estimation (`calculate_mtcn.py`)
Computes the mitochondrial copy number (mtCN) per sample using the mean mitochondrial coverage and the median nuclear coverage. Results are written as a `mtCN_summary.txt` file, which can be used to filter low-quality samples.

### Hail matrix table generation (`annotate_coverage.py`)
This script parses single-sample VCFs together with per-base coverage data to build a Hail MatrixTable (MT). The MT is a unified object that efficiently stores the variants × samples matrix data, including coverage and metadata. This pipeline leverages the pre-calculated per-base coverage to fill in missing information, such as the coverage of a homoplasmic reference site in samples that do not contain a variant at that specific locus.

### Merging VCF files(`combine_vcfs.py`)
This utility performs the merging of multiple VCF files. It takes the per-base coverage of each sample (in the form of a Hail MatrixTable) and other filtering parameters as input. The final output files include a Hail MatrixTable and a compressed gVCF file.

### Functional annotation 
This stage enriches and filters variants using VEP (Variant Effect Predictor) and a comprehensive set of commonly used mtDNA resources.

The functional annotation supports two workflow modes to accommodate different data processing requirements:

1. Combined/Multi-sample Mode 

`add_annotation_combine.py`: Serves as the primary execution script. It first runs VEP on the merged multi-sample gVCF to perform basic gene and consequence annotation. Subsequently, it calls `process_variants_combine.py` to apply deeper-level annotation and filtering.

`process_variants_combine.py`: Reads the VEP annotation output, integrates information from multiple external databases, and performs stringent variant filtering.

2. Single-sample Mode

`add_annotation_single.py`: Responsible for processing single VCF files (typically run in a batch). It runs VEP and then calls process_variants_single.py.

`process_variants_single.py`: Focuses on annotating and filtering variants for a single sample, producing the final single-sample variant list file.

Both process_variants scripts contain the core annotation and filtering logic, integrating data from the following resources:

- ClinVar (germline classification)

- gnomAD mitochondrial frequencies

- HelixMTdb (het/hom AFs and max heteroplasmy)

- MITOMAP (polymorphisms & disease records)

- MitoTIP (tRNA pathogenicity prediction)

- APOGEE / MitImpact (protein impact scoring)

- HmtVar (pathogenicity meta-predictions)

- PhyloTree-based haplogroup variant sets for haplogroup consistency checks
