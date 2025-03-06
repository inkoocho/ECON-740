##
##    Written by Jung Jae Kim on January 7, 2024. Parallel execution code of duopoly_nash_julia_function.jl
##

using Distributed

num_cores = 10;
addprocs(num_cores; exeflags="--project");

@everywhere begin
  using Random, Plots, Distributions, Dates, SharedArrays, LinearAlgebra  
  include("duopoly_nash_julia_function.jl");
  N = 10000;
  lambda = 0.0025;
  sigmawn = 0.0065;
end

simul = 120;
result = SharedArray{Float64}(simul,4,N);

@time @sync @distributed for i in 1:simul  
    res = duopoly_nash_julia_function(N,lambda,sigmawn)
    result[i,:,:] = res
end

using JLD2
@save "duopoly_result.jld"