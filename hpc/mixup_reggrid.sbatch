#!/bin/bash
# Sensitivity check: the split grid REGISTERED in amendment 5 (quantiles
# 0.30-0.70) rather than the widened 0.15-0.85 grid actually used.  Writes
# results/28_mixup_<name>_reg.txt; nothing else changes.
#SBATCH -J ecs-reggrid
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -o joblogs/reggrid_%A_%a.out
#SBATCH --array=0-8
set -e
mkdir -p joblogs
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
export MIXUP_FRACS=0.30,0.40,0.50,0.60,0.70
NAMES=(D1 D2 D3 D4 D5 D6 D7 D8 ELT)
cd experiments
python3 28_mixup_datasets.py "${NAMES[$SLURM_ARRAY_TASK_ID]}"
