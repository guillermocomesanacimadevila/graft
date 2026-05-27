#!/usr/bin/env bash
#SBATCH --job-name=graft
#SBATCH --output=logs/graft_%j.out
#SBATCH --error=logs/graft_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail 

CPUS="${1}"
TIME_HRS="${2}"
MAX_FORKS="${3}"
MEMORY_GB="${4}"
CONTAINER="${5}"
PARAMS_FILE="${6}"

cd "${SLURM_SUBMIT_DIR}"

# load required stuff
mkdir -p logs work env
module purge
module load nextflow
command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not found"; exit 1; }

nextflow run . \
  -profile singularity \
  -params-file "${PARAMS_FILE}" \
  --container "${CONTAINER}" \
  --cpus "${CPUS}" \
  --time_hrs "${TIME_HRS}" \
  --memory "${MEMORY_GB} GB" \
  --max_forks "${MAX_FORKS}" 
