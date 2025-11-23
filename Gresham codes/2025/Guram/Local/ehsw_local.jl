using Distributed
using CSV, DataFrames

num_cores = 10;
addprocs(num_cores; exeflags="--project");

@everywhere begin
  using Random, Plots, Distributions, Dates, SharedArrays
  include("ehsw_function_local.jl");
  Titer = 6000;
  alpha = 0.96;
  rho = 0.95;
  sigmap = 1.0;
  sigmaf = 1.0;
  sigmab = 0.00001;    
end

simul = 120;
#result = SharedArray{Float64}(simul,Titer,4);

results = @time @sync @distributed (vcat) for i in 1:simul
  ehsw_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)
end

# Reshape each simulation result

reshaped_results = [reshape(results[(i-1)*Titer+1:i*Titer, :], 1, Titer, 4) for i in 1:simul]

# Concatenate the results along the first dimension to create a 3D array
final_result = vcat(reshaped_results...)

#@time @sync @distributed for i in 1:simul  
#    res = gresham4_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)
#    result[i,:,:] .= res
#end

#=
for i in 1:simul
    plot1 = plot(result[i,:,1],label="Price")
    plot2 = plot(result[i,:,2],ylim=(0,1),label="Pi")
    plot3 = plot(result[i,:,3],label="Beta0")
    plot4 = plot(result[i,:,4],label="Beta1")
    p = plot(plot1, plot2, plot3, plot4, layout = (2, 2))
    display(p)
end
=#

#simul, Titer, variable =size(result) 
b_range=range(0,1.1,length=41)
histogram(final_result[:,Titer,2], label="terminal prob", bins =b_range)




# Flatten final_result (simul, Titer, 4) into a table with columns: sim, time, price, pi, beta0, beta1
nsim, ntime, nvar = size(final_result)
sim_ids = repeat(1:nsim, inner=ntime)
times = repeat(1:ntime, outer=nsim)
vals = reshape(final_result, nsim*ntime, nvar)

df = DataFrame(sim = sim_ids,
         time = times,
         price = vals[:,1],
         pi = vals[:,2],
         beta0 = vals[:,3],
         beta1 = vals[:,4])

# Write CSV to same directory as this script
outdir = dirname(@__FILE__)
outpath = joinpath(outdir, "ehsw_path.csv")
CSV.write(outpath, df)
println("Wrote ehsw_path to ", outpath)