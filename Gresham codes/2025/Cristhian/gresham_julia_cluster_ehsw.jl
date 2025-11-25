using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere begin
    include("gresham_julia_function_ehsw.jl");
    Titer = 20000;
    alpha = 0.96;
    rho = 0.95;
    sigmap = 1.0;
    sigmaf = 1.0;
    sigmab = 0.00001; 
end

simul = 10000;

results = @time @sync @distributed (vcat) for i in 1:simul
  gresham_julia_function_ehsw(Titer,alpha,rho,sigmap,sigmaf,sigmab)
end

using JLD2
@save "gresham_result.jld"
