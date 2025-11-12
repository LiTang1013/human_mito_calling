#!/usr/bin/env bash

#SBATCH -J env_setup
#SBATCH -o env_setup_%j.out
#SBATCH -e env_setup_%j.err
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00

# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error.
set -eo pipefail

echo "========== Job Started =========="
date

# ---------------- Locate Conda Base ----------------
# Temporarily disable 'nounset' (-u) before and after sourcing conda.sh
# to prevent issues with unbound variables within the environment setup scripts.

# Prioritize the user-provided path; fallback to 'conda info --base'
CONDA_BASE="/path/to/conda_base" # Change this path if needed.

if [[ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # Use 'conda info --base' only if 'conda' command is found
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  fi
fi

if [[ -z "${CONDA_BASE}" || ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  echo "[ERROR] conda.sh not found. Set CONDA_BASE correctly or ensure 'conda' is on PATH." >&2
  exit 1
fi

# Source Conda initialization script
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

# ---------------- Ensure Mamba is installed in base ----------------
if ! command -v mamba >/dev/null 2>&1; then
  echo "[INFO] Installing mamba into base environment..."
  set +u
  conda activate base
  set -u
  conda install -n base -c conda-forge -y mamba
  conda deactivate
fi

ENV_NAME="bio_pipeline_env"

# ---------------- Create or Update Environment ----------------
# Check if environment already exists
if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[INFO] Creating environment '${ENV_NAME}'..."
  # Package list:
  # - openjdk=17 for GATK/Picard/Hail
  # - ensembl-vep=113.0
  # - Perl libs (perl-list-moreutils, perl-json, etc.) for VEP dependencies
  mamba create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
    "python=3.11" \
    "openjdk=17" \
    "gatk4=4.5.0.0" \
    "picard-slim=3.1.1" \
    "bwa=0.7.18" \
    "samtools=1.21" \
    "r-base=4.3" \
    "perl=5.32.*" \
    "perl-list-moreutils=0.430.*" \
    "perl-json" \
    "perl-archive-zip" \
    "perl-dbi" \
    "perl-dbd-mysql" \
    "ensembl-vep=113.0"
else
  echo "[INFO] Environment '${ENV_NAME}' already exists. Skipping creation."
fi

# Activate the environment
set +u
conda activate "${ENV_NAME}"
set -u

# ---------------- Install Python Libraries ----------------
python -m pip install --upgrade pip
# Install Hail and matching dependencies. Letting Hail pull its own PySpark/py4j
# minimizes manual dependency conflicts.
python -m pip install \
  "hail==0.2.133" \
  "pandas>=2.2,<3" \
  "numpy>=1.26,<3"

# ---------------- Sanity Checks ----------------
echo "---------- Sanity Checks ----------"

# Java Checks
which java || true
java -version || true

# GATK Check
echo "gatk version:"
gatk --version || true

# Picard Check (Suppressing full output)
echo "picard help check..."
picard -h >/dev/null || true

# VEP Check
echo "vep version:"
vep --help | head -n 5 || true

# Perl Dependency Check (List::MoreUtils)
perl -e 'use List::MoreUtils qw(any all); print "Perl List::MoreUtils OK\n";' || true

# Python and Hail Checks
python -V
python - <<'PY'
try:
    import hail as hl
    import pyspark
    import py4j
    print("[OK] Hail version:", hl.__version__)
    print("[INFO] PySpark version:", pyspark.__version__)
    print("[INFO] Py4J version:", py4j.__version__)
except Exception as e:
    print("[WARN] Hail import failed:", e)
PY

# Export for Spark/Hail (Optional but good practice)
export PYSPARK_PYTHON="$(command -v python)"
echo "PYSPARK_PYTHON=${PYSPARK_PYTHON}"

echo "========== Job Finished =========="
date