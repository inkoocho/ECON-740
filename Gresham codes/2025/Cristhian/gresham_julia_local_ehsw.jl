##
##    Written by Jung Jae Kim on December 29, 2023. Parallel execution code of gresham4_julia_function.jl
## Modified by Cristhian Rosales-Castillo to estimate the Gresham model using Cobweb model 
## Last modified: November 24, 2025
## Class: Learning and AI in Economics
## Professor: In-Koo Cho

using Distributed

num_cores = 6;
addprocs(num_cores; exeflags="--project");

@everywhere begin
  using Random, Plots, Distributions, Dates, SharedArrays
  include("gresham_julia_function_ehsw.jl");
  Titer = 2000;
  alpha = 0.96;
  rho = 0.95;
  sigmap = 1.0;
  sigmaf = 1.0;
  sigmab = 0.00001;    
end

    simul = 50;
    #result = SharedArray{Float64}(simul,Titer,4);

    results = @time @sync @distributed (vcat) for i in 1:simul
    gresham_julia_function_ehsw(Titer,alpha,rho,sigmap,sigmaf,sigmab)
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

#using JLD2
#@save "gresham_result.jld"
