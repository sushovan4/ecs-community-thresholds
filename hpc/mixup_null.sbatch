#!/bin/bash
#SBATCH -J ecs-mixnull
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o mixnull_%j.out
set -e
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NPROC=$SLURM_CPUS_PER_TASK
export MIXUP_PATH="$SLURM_SUBMIT_DIR/vendor"
cd experiments
python3 30_mixup_largeP_null.py
