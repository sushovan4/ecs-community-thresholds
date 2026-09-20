#!/bin/bash
#SBATCH -J ecs-temp-occ
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o temporal_occ_%A_%a.out
#SBATCH --array=0-2
set -e
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SERIES=(T2 T3 T4)
cd experiments
python3 22_temporal.py run "${SERIES[$SLURM_ARRAY_TASK_ID]}"
