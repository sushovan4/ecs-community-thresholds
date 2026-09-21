#!/bin/bash
# Amendment 8 inputs: seeded TITAN with named indicators, one dataset per task.
# R comes from a user-space conda env because the site module is a broken
# container (see hpc/README.md).
#SBATCH -J ecs-titan
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 48:00:00
#SBATCH -o joblogs/titan_%A_%a.out
#SBATCH --array=0-2
set -e
cd "$SLURM_SUBMIT_DIR"
mkdir -p joblogs
module load anaconda/2023.03
source activate "$HOME/renv"
NAMES=(D5 D8 D1)
cd experiments
Rscript 35_titan_sppmax.R "${NAMES[$SLURM_ARRAY_TASK_ID]}"
