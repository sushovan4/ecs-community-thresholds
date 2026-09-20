#!/bin/bash
#SBATCH -J ecs-temporal
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o temporal_%A_%a.out
#SBATCH --array=0-9
#
# Temporal precedence, primary test per series (preregistration/temporal-precedence.md).
#   cd ~/HYPERDESIGN && git pull && sbatch hpc/temporal.sbatch
set -e
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SERIES=(E2 E6 E3 E4 E8 C1 C2 C3 C4 N1)
cd experiments
python3 22_temporal.py run "${SERIES[$SLURM_ARRAY_TASK_ID]}"
