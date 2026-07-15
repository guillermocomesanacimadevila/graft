#!/usr/bin/env bash
set -euo pipefail

# this will be executed from the main Jupyter notebook
# * log into SLURM cluster
# * load Apptainer and Nextflow
# * clone / update graft
# * check whether the Apptainer container is already there
# * if not, pull it from GHCR


REPO_URL="https://github.com/guillermocomesanacimadevila/graft.git"
REPO_DIR="$HOME/graft"
SIF_DIR="$REPO_DIR/env"
SIF_FILE="$SIF_DIR/graft_1.0.0.sif"
IMAGE_URI="docker://ghcr.io/guillermocomesanacimadevila/graft:1.0.0"

# load HPC modules
if command -v module >/dev/null 2>&1; then
    module load apptainer
    module load nextflow
fi


# check required software
command -v apptainer >/dev/null 2>&1 || {
    echo "ERROR: apptainer not found"
    exit 1
}

command -v nextflow >/dev/null 2>&1 || {
    echo "ERROR: nextflow not found"
    exit 1
}

command -v git >/dev/null 2>&1 || {
    echo "ERROR: git not found"
    exit 1
}


# clone / update graft
if [ ! -d "$REPO_DIR/.git" ]; then
    echo "[TRACKING] graft repo not found. Cloning..."
    git clone "$REPO_URL" "$REPO_DIR"
else
    echo "[TRACKING] graft repo already exists. Updating..."
    git -C "$REPO_DIR" pull --ff-only
fi


# create container directory
mkdir -p "$SIF_DIR"


# pull container if not already present
if [ -f "$SIF_FILE" ]; then
    echo "[TRACKING] Apptainer image already exists: $SIF_FILE"
else
    echo "[TRACKING] Pulling Apptainer image from GHCR..."
    apptainer pull "$SIF_FILE" "$IMAGE_URI"
fi


# confirm container works
echo "[TRACKING] Checking Apptainer image..."
apptainer exec "$SIF_FILE" python3 --version
apptainer exec "$SIF_FILE" Rscript --version


echo "[DONE] graft HPC environment ready."
echo "[INFO] Repo: $REPO_DIR"
echo "[INFO] SIF:  $SIF_FILE"