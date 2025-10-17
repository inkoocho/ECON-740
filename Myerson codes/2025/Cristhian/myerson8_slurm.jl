# ********************************************
# Model simulation of Cho and Libgober (2025)
# Class: Learning and AI in Economics, Emory University
# Professor: In-Koo Cho
# Date: October 16, 2025
# Name: Cristhian Rosales-Castillo
# ********************************************

using Distributed, SlurmClusterManager

addprocs(SlurmManager())

@everywhere begin 
    using Distributions, Roots, LinearAlgebra, Random, Statistics, CSV, DataFrames
    include("myerson8_function.jl")
    include("truncated_normal_pdf.jl")
    include("truncated_normal_cdf.jl")
    include("hazard_rate.jl")
    include("truncated_normal_sample.jl")
end

mu = 10
sigmamin = 11
sigmaincre = 0.001
nsigma = 3000

sigmavec = [sigmamin + sigmaincre * (j - 1) for j in 1:nsigma]
N = 100
a = 0.0001
epsilon = 0.75
T = 300000
CRsample = 10

sims = @time pmap(1:nsigma) do isigma
    myerson8_function(mu, sigmavec[isigma], N, a, epsilon, T, CRsample)
end

CLferror_vec = [sims[i][1] for i in 1:nsigma]
CRferror_mat = hcat([sims[i][2] for i in 1:nsigma]...)

CSV.write("CLferror_vec.csv", DataFrame(CLferror = CLferror_vec))
CSV.write("CRferror_mat.csv", DataFrame(CRferror_mat, :auto))
