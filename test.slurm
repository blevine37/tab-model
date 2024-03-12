#!/bin/bash
#
#SBATCH --job-name=exa
#SBATCH --output=res.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p long-28core
#SBATCH --mail-user=fangchun.liang@stonybrook.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

gfortran exact.f -o run.exe
./run.exe 
