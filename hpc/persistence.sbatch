#!/bin/bash
#SBATCH -J ecs-persistence
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o persistence_%j.out
#
# "Why not persistence?": persistence landscapes and total persistence vs chi
# through the Table 2 estimator (experiments/20_persistence_swap.py).
# Small request (8 cores, 8 h) so the scheduler can backfill it quickly.
#   cd ~/HYPERDESIGN && git pull && sbatch hpc/persistence.sbatch
set -e
cd "$SLURM_SUBMIT_DIR"
module load python3/3.13.3      # the venv links against this build
source venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NPROC=$SLURM_CPUS_PER_TASK
cd experiments
python3 20_persistence_swap.py
