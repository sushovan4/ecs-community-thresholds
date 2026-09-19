#!/bin/bash
#SBATCH -J hyperdesign-sim
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o hyperdesign_%j.out
#
# Publication-strength simulation runs for paper sections 3.3-3.5.
# Submit from ~/HYPERDESIGN after the one-time venv setup:  sbatch hpc/pegasus_full.sbatch
set -e
cd "$SLURM_SUBMIT_DIR"
source venv/bin/activate
cd experiments
export NPERM=499 NBOOT=200
python3 ../hpc/driver_full.py
