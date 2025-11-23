

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere begin
    include("ehsw_function_local.jl");
    Titer = 20000;
    alpha = 0.96;
    rho = 0.95;
    sigmap = 1.0;
    sigmaf = 1.0;
    sigmab = 0.00001; 
end

simul = 10000;

results = @time @sync @distributed (vcat) for i in 1:simul
  ehsw_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)
end


using CSV, DataFrames



df = DataFrame(sim = Int[], time = Int[], price = Float64[], pi = Float64[], beta0 = Float64[], beta1 = Float64[])
nsim = simul
ntime = Titer
for sim_id in 1:nsim
    for time in 1:ntime
        idx = (sim_id - 1) * ntime + time
        push!(df, (sim_id, time, results[idx, 1], results[idx, 2], results[idx, 3], results[idx, 4]))
    end
end
outpath = "ehsw_cpu_results.csv"

CSV.write(outpath, df)
