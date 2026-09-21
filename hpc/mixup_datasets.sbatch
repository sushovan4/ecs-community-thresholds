#!/bin/bash
#SBATCH -J ecs-mixdata
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -o mixdata_%A_%a.out
#SBATCH --array=0-6
set -e
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
NAMES=(ELT T2 T3 T4 E2 E3 E4)
cd experiments
python3 28_mixup_datasets.py "${NAMES[$SLURM_ARRAY_TASK_ID]}"
