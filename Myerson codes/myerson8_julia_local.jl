using Distributed

num_cores = 20;
addprocs(num_cores; exeflags="--project");

@everywhere begin
    using Distributions, Random, Statistics, LinearAlgebra, Roots
    include("myerson8_julia_function.jl");
    mu = 10              
    sigmamin = 11
    sigmaincre = 0.001    
    Nbuyer = 100
    a = 0.0001
    T = 300000
    epsilon = 0.75
    CRsample = 5
end

jsigma = 200  
results = @time @sync @distributed (vcat) for j in 1:jsigma
  sigma = sigmamin + sigmaincre * (j - 1)
  myerson8_julia_function(mu, sigma, Nbuyer, a, epsilon, T, CRsample)
end

CLferror = results[:,1]
CRferror = results[:,2:CRsample+1]
