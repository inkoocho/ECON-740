#!/bin/tcsh

#SBATCH -p debug
#SBATCH --time=0
#SBATCH --mem=MaxMemPerNode
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=20
#SBATCH --exclude=inkoocho
#SBATCH --job-name=test
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

julia myerson8_julia_cluster.jl


