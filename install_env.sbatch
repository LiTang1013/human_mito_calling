#!/usr/bin/env bash
#SBATCH -J env_setup
#SBATCH -o env_setup_%j.out
#SBATCH -e env_setup_%j.err
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00

# 严格模式，但在激活 conda 前后会临时关闭 -u
set -euo pipefail

echo "========== Job Started =========="
date

# --------- 定位 Conda base ----------
# 优先使用你提供的路径；若没改或无效，自动 fallback 到 `conda info --base`
CONDA_BASE="/apps/software/2022b/software/miniconda/24.11.3"  # 如需更改，改这里
if [[ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  fi
fi

if [[ -z "${CONDA_BASE}" || ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  echo "[ERROR] conda.sh not found. Set CONDA_BASE correctly or ensure 'conda' is on PATH." >&2
  exit 1
fi

# 规避 activate.d 与 set -u 冲突
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

# 确保 base 里有 mamba
if ! command -v mamba >/dev/null 2>&1; then
  echo "[INFO] Installing mamba into base..."
  set +u
  conda activate base
  set -u
  conda install -n base -c conda-forge -y mamba
  conda deactivate
fi

ENV_NAME="bio_pipeline_env"

# 若环境不存在则创建（存在则跳过创建）
if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[INFO] Creating env '${ENV_NAME}' ..."
  # 说明：
  # - openjdk=17 供 GATK/Picard/Hail
  # - ensembl-vep=113.0
  # - perl-list-moreutils=0.430.* 解决 VEP Perl 依赖
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
fi

# 激活环境（同样临时关闭 -u 以防 activate.d 报 unbound variable）
set +u
conda activate "${ENV_NAME}"
set -u

# --------- Python libs ----------
python -m pip install --upgrade pip
# 让 hail 自动拉匹配的 pyspark/py4j，避免手工指定引发的依赖冲突
python -m pip install \
  "hail==0.2.133" \
  "pandas>=2.2,<3" \
  "numpy>=1.26,<3"

# --------- Sanity checks ----------
python - <<'PY'
try:
    import hail as hl, pyspark, py4j
    print("hail   =", hl.__version__)
    print("pyspark=", pyspark.__version__)
    print("py4j   =", py4j.__version__)
except Exception as e:
    print("[WARN] Python libs check failed:", e)
PY

echo "---------- Sanity Checks ----------"
which java || true
java -version || true

echo "gatk version:"
gatk --version || true

echo "picard help:"
picard -h >/dev/null || true

echo "vep version:"
vep --help | head -n 5 || true
# 额外：验证 Perl 依赖（List::MoreUtils）
perl -e 'use List::MoreUtils qw(any all); print "Perl List::MoreUtils OK\n";' || true

echo "python:"
python -V
python - <<'PY'
try:
    import hail as hl
    print("[OK] Hail version:", hl.__version__)
except Exception as e:
    print("[WARN] Hail import failed:", e)
PY

# 为 Spark 指定解释器（可选）
export PYSPARK_PYTHON="$(command -v python)"
echo "PYSPARK_PYTHON=${PYSPARK_PYTHON}"

echo "========== Job Finished =========="
date
