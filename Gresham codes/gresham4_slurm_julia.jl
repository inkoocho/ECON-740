using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere begin
    include("gresham4_julia_function.jl");
    Titer = 6000;
    alpha = 0.96;
    rho = 0.95;
    sigmap = 1.0;
    sigmaf = 1.0;
    sigmab = 0.00001; 
end

simul = 100;

results = @time @sync @distributed (vcat) for i in 1:simul
  gresham4_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)
end

using JLD2
@save "gresham_result.jld"