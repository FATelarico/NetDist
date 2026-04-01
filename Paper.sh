#!/bin/bash
#SBATCH --job-name=JIBE
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --output=JIBE.out

# Absolute paths
SIF="$HOME/Apptainer/graph-tool.sif"
SCRIPT="$HOME/JIBE/Paper.py"

# Go to a working directory where output should appear
cd "$HOME/JIBE"


# Run inside container
srun apptainer exec "$SIF" python3 "$SCRIPT"
