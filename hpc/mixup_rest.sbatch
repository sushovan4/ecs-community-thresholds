#!/bin/bash
# Amendment 5 second pass, the items not yet run: the remaining registered
# temporal series and the three spatial datasets D6-D8.
#SBATCH -J ecs-mixrest
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 16:00:00
#SBATCH -o joblogs/mixrest_%A_%a.out
#SBATCH --array=0-9
set -e
cd "$SLURM_SUBMIT_DIR"
mkdir -p joblogs
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
NAMES=(D6 D7 D8 E6 E8 C1 C2 C3 C4 N1)
cd experiments
python3 28_mixup_datasets.py "${NAMES[$SLURM_ARRAY_TASK_ID]}"
