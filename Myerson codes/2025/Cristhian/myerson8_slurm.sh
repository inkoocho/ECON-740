#!/bin/bash

#SBATCH -p econ
#SBATCH --time=0-3:00:00
#SBATCH --mem=64G
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=myerson
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

module load julia 

julia myerson8_slurm.jl