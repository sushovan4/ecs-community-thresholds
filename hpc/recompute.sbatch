#!/bin/bash
# Recomputation under GUDHI precision="safe" (the 2026-09-20 correction).
# Every interaction-test number in the paper, widened and registered grids.
#SBATCH -J ecs-recomp
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH -t 48:00:00
#SBATCH -o joblogs/recomp_%A_%a.out
#SBATCH --array=0-21
set -e
cd "$SLURM_SUBMIT_DIR"
mkdir -p joblogs
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
NAMES=(D1 D2 D3 D4 D5 D6 D7 D8 ELT E2 E3 E4 E6 E8 C1 C2 C3 C4 N1 T2 T3 T4)
D="${NAMES[$SLURM_ARRAY_TASK_ID]}"
cd experiments
python3 28_mixup_datasets.py "$D"
case "$D" in
  D*|ELT) MIXUP_FRACS=0.30,0.40,0.50,0.60,0.70 python3 28_mixup_datasets.py "$D" ;;
esac
