# Econ 777 

Matlab and Julia code repositories for ECON 777 class

## 1. Gresham codes 
The Matlab and Julia codes are based on [Cho and Kasa (2017, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20160665).

- `gresham4_matlab_forparallel.m`
    - Matlab function script that generates sample paths of $\pi_t$, $\beta_t(0)$, $\beta_t(1)$, $p_t$
- `gresham4_julia_function.jl`
    - Julia function script translated from `gresham4_matlab_forparallel.m`.
    - It is slightly modified to have parameters as inputs
- `gresham4_julia_local.jl`
    - Parallel execution code of `gresham4_julia_function.jl` on your **local** computer.
    - Adjust 'num_cores' depending on your CPU.
    - This code also generates some plots and histograms from the simulation results.
- `gresham4_julia_cluster.jl`
    - Parallel execustion code of `gresham4_julia_function.jl` on the **cluster**.
    - Can be executed using the shell file `gresham4_slurm.sh`.

## 2. Duopoly codes
The Matlab and Julia codes are based on [this paper](https://github.com/jay9209/ECON-777/blob/main/Duopoly%20codes/duopoly4.pdf)

## 3. Compete codes
The Matlab and Julia codes are based on [this paper](https://github.com/jay9209/ECON-777/blob/main/Compete%20codes/competing10.pdf)