#!/usr/bin/env python3
import subprocess
import datetime
import os 
import argparse 
from pathlib import Path

# steps
# parse SBATCH params
# > CPUs
# > t
# > GPUs
# container build (docker / singularity)
# assets/ dir in HPC
# slurm id 
# slurm or local 

def perform_checks(container_build: str, slurm_id: str, assets_dir: str, platform: str, cpus: int, time: int, maxForks: int, memory: int) -> dict:
    
    """
    Validate user inputs before building the nf-core/graft command.
    Returns cleaned values if everything looks good.
    """

    allowed_containers = {"docker", "singularity"}
    allowed_platforms = {"local", "slurm"}

    container_build = container_build.strip().lower()
    slurm_id = slurm_id.strip()
    platform = platform.strip().lower()
    assets_dir = Path(assets_dir).expanduser().resolve()

    if not container_build:
        raise ValueError("Empty container build. Choose either: docker or singularity.")

    if container_build not in allowed_containers:
        raise ValueError(
            f"Invalid container build: {container_build}. "
            "Choose either: docker or singularity."
        )

    if not platform:
        raise ValueError("Empty platform. Choose either: local or slurm.")

    if platform not in allowed_platforms:
        raise ValueError(
            f"Invalid platform: {platform}. Choose either: local or slurm."
        )

    if not assets_dir.exists():
        raise FileNotFoundError(f"Assets directory does not exist: {assets_dir}")

    if not assets_dir.is_dir():
        raise NotADirectoryError(f"Assets path is not a directory: {assets_dir}")

    required_assets = ["params.stage1.yaml", "gwas.tsv", "ldsc_pairs.tsv"]
    missing = []
    for file in required_assets:
        if not (assets_dir / file).exists():
            missing.append(file)

    if missing:
        raise FileNotFoundError(
            "Missing required files in assets directory: "
            + ", ".join(missing)
        )

    if platform == "slurm" and not slurm_id:
        raise ValueError("Come on big man, include a proper HPC / SLURM account ID.")

    if cpus <= 0:
        raise ValueError("CPUs must be greater than 0.")

    if time <= 0:
        raise ValueError("Time must be greater than 0.")

    if maxForks <= 0:
        raise ValueError("maxForks must be greater than 0.")

    if memory <= 0:
        raise ValueError("Memory must be greater than 0.")

    if maxForks > cpus:
        raise ValueError(
            f"maxForks ({maxForks}) cannot be greater than CPUs ({cpus})."
        )

    return {
        "container_build": container_build,
        "slurm_id": slurm_id,
        "assets_dir": assets_dir,
        "platform": platform,
        "cpus": cpus,
        "time": time,
        "maxForks": maxForks,
        "memory": memory,
    }

# nextflow run . \
#   -profile docker \
#   -c conf/local/nextflow.config \
#   -params-file assets/params.stage1.yaml \
#   --input assets/gwas.tsv \
#   --pairs assets/ldsc_pairs.tsv \
#   -process.maxForks 1

# nextflow.config == slurm - platform
# - platform slurm or local and based on that we select the nextflow.config dir 
# maxForks


# We need to pull singularity container in HPC
# and do an if statement to check whether its already been built
# if not then build 

def prep_slurm_log_and_query(
    container_build: str,
    slurm_id: str,
    cpus: int,
    time: int,
    assets_dir: str,
    platform: str,
    maxForks: int,
    memory: int
):
    """
    Validate inputs locally, then SSH into Falcon, ensure the Apptainer
    container exists, and submit the Nextflow SLURM job.
    """

    config = perform_checks(
        container_build=container_build,
        slurm_id=slurm_id,
        assets_dir=assets_dir,
        time=time,
        cpus=cpus,
        maxForks=maxForks,
        platform=platform,
        memory=memory
    )

    slurm_id           = config["slurm_id"]
    memory             = config["memory"]
    cpus               = config["cpus"]
    time               = config["time"]
    maxForks           = config["maxForks"]
    remote_project_dir = "~/graft"
    sif                = "env/graft_1.0.0.sif"
    image              = "docker://ghcr.io/guillermocomesanacimadevila/graft:1.0.0"
    params_file        = "assets/params.stage1.yaml"

    remote_cmd = f"""
set -euo pipefail

cd {remote_project_dir}
mkdir -p env logs work
command -v apptainer >/dev/null 2>&1 || {{ echo "ERROR: apptainer not found"; exit 1; }}
if [[ ! -f "{sif}" ]]; then
    echo "Container not found: {sif}"
    echo "Pulling: {image}"
    apptainer pull "{sif}" "{image}"
else
    echo "Container already exists: {sif}"
fi

sbatch \\
  --cpus-per-task={cpus} \\
  --mem={memory}G \\
  --time={time}:00:00 \\
  bin/slurm.sh \\
  {cpus} \\
  {time} \\
  {maxForks} \\
  {memory} \\
  "{sif}" \\
  "{params_file}"
"""
    # Log into falcon cluster 
    subprocess.run(["ssh", f"{slurm_id}@falconlogin.cf.ac.uk", remote_cmd], check=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--container_build", type=str, required=True)
    parser.add_argument("--slurm_id", type=str, required=False, default="")
    parser.add_argument("--assets_dir", type=str, required=True)
    parser.add_argument("--platform", type=str, required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--time", type=int, required=True)
    parser.add_argument("--maxForks", type=int, required=True)
    parser.add_argument("--memory", type=int, required=True)
    args = parser.parse_args()
    if args.platform.lower() == "slurm":
        prep_slurm_log_and_query(
            container_build=args.container_build,
            slurm_id=args.slurm_id,
            cpus=args.cpus,
            time=args.time,
            assets_dir=args.assets_dir,
            platform=args.platform,
            maxForks=args.maxForks,
            memory=args.memory,
        )

    else:
        raise NotImplementedError(
            "Local execution pathway not implemented yet."
        )

if __name__ == "__main__":
    main()
