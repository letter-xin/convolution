#!/bin/bash

#SBATCH -J xfx
#SBATCH -p cpu
#SBATCH -n 1
##SBATCH --ntasks-per-node=1
#SBATCH -o output/system/%j.out
#SBATCH -e output/system/%j.err
#SBATCH --time=48:00:00

ulimit -l unlimited
ulimit -s unlimited

./ILS