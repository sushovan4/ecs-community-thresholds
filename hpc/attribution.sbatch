#!/bin/bash
# Amendment 8 inputs: seeded TITAN with named indicators, and the curvature
# attribution, for every spatial dataset.
#SBATCH -J ecs-attrib
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o joblogs/attrib_%A_%a.out
#SBATCH --array=0-7
set -e
cd "$SLURM_SUBMIT_DIR"
mkdir -p joblogs
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
NAMES=(D7 D6 D2 D3 D4 D5 D1 D8)
D="${NAMES[$SLURM_ARRAY_TASK_ID]}"
cd experiments
python3 36_attribution_multi.py "$D"
# TITAN runs locally: R on Pegasus is containerized and fails to mount
