#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --time=0-1:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus-per-task=nvidia_a30_2g.12gb:1
#SBATCH -A econ740s1
#SBATCH -e slurm.%j.err

export SLURM_EXPORT_ENV=ALL
export PATH=/usr/bin

# srun
PATH=/opt/slurm/default/bin:"${PATH}"

# julia
PATH=/scion/julia/1.10.9-bin/bin:"${PATH}"

# cuda
export CUDA_HOME=/scion/cuda/12.6

echo '=BEG=' "$(date)"

# SlurmClusterManager invokes srun
srun julia ehsw_gpu.jl

echo '=END=' "$(date)"

sleep 30