#!/bin/tcsh

#SBATCH -p debug
#SBATCH --time=0
#SBATCH --mem=MaxMemPerNode
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --exclude=inkoocho
#SBATCH --job-name=test
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

octave -q --eval 'myerson8_matlab.parallel' &
