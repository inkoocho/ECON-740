#!/bin/tcsh

#SBATCH -p debug
#SBATCH --time=0
#SBATCH --mem=MaxMemPerNode
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=100
#SBATCH --exclude=inkoocho
#SBATCH --job-name=test
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

julia compete_julia_pihat_cluster.jl
