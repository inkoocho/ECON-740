#############################
# Replication of Cho-Libgober
#############################

using Distributed

addprocs()

@everywhere using Distributions, Random, Statistics, LinearAlgebra, DataFrames, CSV, Roots



# Truncated-at-a CDF: F_trunc(x | X ≥ a)
@everywhere function truncated_normal_a_cdf(x::Float64, mu::Float64, sigma::Float64, a::Float64)
    if x < a
        return 0.0
    else
        α = (a - mu) / sigma
        ξ = (x - mu) / sigma
        Φα = cdf(Normal(0,1), α)
        Φξ = cdf(Normal(0,1), ξ)
        return (Φξ - Φα) / (1.0 - Φα)
    end
end

# Inverse-CDF sampler from Normal(μ,σ) truncated below at a
@everywhere function truncated_normal_a_sample(mu::Float64, sigma::Float64, a::Float64, n::Int)
    α = (a - mu) / sigma
    Φα = cdf(Normal(0,1), α)
    u  = rand(n)                          # U ~ Unif(0,1)
    ξcdf = Φα .+ u .* (1.0 - Φα)          # map to [Φ(α), 1]
    ξ = quantile.(Ref(Normal(0,1)), ξcdf) # standard normal inverse CDF
    return mu .+ sigma .* ξ
end

# Correct p*: solve g(p) = (1 - F(p)) - p * f(p) = 0 with ORIGINAL Normal(μ,σ)
@everywhere function pstar_true(mu::Float64, sigma::Float64)
    d = Normal(mu, sigma)
    g(p) = (1.0 - cdf(d, p)) - p * pdf(d, p)
    # Robust bracketing around μ; the sign change is near μ 
    lo = mu
    hi = mu + 50.0                         
    
    vlo = g(lo); vhi = g(hi)
    expand = 0
    while vlo * vhi > 0 && expand < 5
        hi += 50.0
        vhi = g(hi)
        expand += 1
    end
    
    find_zero(g, (lo, hi), Bisection(); xatol=1e-10, rtol=1e-10, maxevals=200_000)
end



@everywhere function simulate_sigma(j::Int, mu::Float64, sigmamin::Float64, sigmaincre::Float64,
                                    Nbuyer::Int, a_learn::Float64, epsilon::Float64,
                                    T::Int, CRsample::Int; save_traj::Bool=false)

    sigma = sigmamin + sigmaincre*(j-1)

    # --- True optimum for this (μ,σ) ---
    pstar = try
        pstar_true(mu, sigma)
    catch
        NaN
    end

    # Optimal revenue under truncation at a=μ
    pistar = isfinite(pstar) ? pstar * (1.0 - truncated_normal_a_cdf(pstar, mu, sigma, mu)) : NaN

    # Adaptive learning simulation
    beta = zeros(2, T+1)
    beta[:,1] = [10.0, -5.0]
    e = (rand(Bernoulli(0.5), T) .== 1) .* (2*epsilon) .- epsilon

    eprice_last = NaN
    eprice_traj = save_traj ? Vector{Float64}(undef, T) : Float64[]

    for t in 1:T
        β1, β2 = beta[1,t], beta[2,t]
        eprice = -β1 / (2*β2)                 # belief-optimal price at t
        p      = eprice + e[t]

        if save_traj
            eprice_traj[t] = eprice
        end

        # Demand under truncation at μ
        Fp   = truncated_normal_a_cdf(p, mu, sigma, mu)
        Dagg = rand(Binomial(Nbuyer, 1.0 - Fp))
        q    = Dagg / Nbuyer

        # R_t matrix and update
        R = @views [1.0 eprice; eprice eprice^2 + epsilon^2]
        if q == 0.0
            beta[:, t+1] = beta[:, t] .+ a_learn .* ([-epsilon, 0.0] .- beta[:, t])
        elseif q == 1.0
            beta[:, t+1] = beta[:, t] .+ a_learn .* ([ epsilon, 0.0] .- beta[:, t])
        else
            ϕ   = [1.0, p]
            err = q - dot(ϕ, beta[:, t])
            step = R \ (ϕ * err)
            beta[:, t+1] = beta[:, t] .+ a_learn .* step
        end

        eprice_last = eprice
    end

    # CL forecast error at final period
    CLferror = isfinite(pstar) ? (eprice_last - pstar) : NaN

    # ColeRoughgarden sampling benchmark
    CRferror    = zeros(CRsample)
    estimatedcr = zeros(CRsample)
    for isample in 1:CRsample
        Obvalue = 2 * T * isample
        value   = truncated_normal_a_sample(mu, sigma, mu, Obvalue)
        pcr     = sort(value, rev=true)
        qgrid   = (1:Obvalue) ./ Obvalue
        rev     = qgrid .* pcr
        _, I    = findmax(rev)
        p_hat   = pcr[I]
        estimatedcr[isample] = p_hat
        CRferror[isample]    = isfinite(pstar) ? (pstar - p_hat) : NaN
    end

    return (sigma, pstar, pistar, CLferror, CRferror, estimatedcr, eprice_last, beta[:,end], eprice_traj)
end

# ---------------------------
# Parameters
# ---------------------------
mu         = 10.0
sigmamin   = 11.0
sigmaincre = 0.001
jsigma     = 5000
Nbuyer     = 100
a_learn    = 1e-4
epsilon    = 0.75
T          = 300_000
CRsample   = 5

const SAVE_TRAJ = false   

println("Workers: ", nworkers())


results = pmap(1:jsigma) do j
    simulate_sigma(j, mu, sigmamin, sigmaincre, Nbuyer, a_learn, epsilon, T, CRsample; save_traj=SAVE_TRAJ)
end

# ---------------------------
# CSV
# ---------------------------
cols = ["mu","sigma","pstar_true","pistar_true","CLferror",
        "CRferror1","CRferror2","CRferror3","CRferror4","CRferror5",
        "estcr1","estcr2","estcr3","estcr4","estcr5",
        "eprice_T","beta1_T","beta2_T"]

syms = Symbol.(cols)
rows = NamedTuple[]
for (j, r) in enumerate(results)
    sigma, pstar, pistar, CLferror, CRferror, estcr, eprice_T, beta_T, eprice_traj = r
    push!(rows, NamedTuple{Tuple(syms)}((
        mu, sigma, pstar, pistar, CLferror,
        CRferror[1], CRferror[2], CRferror[3], CRferror[4], CRferror[5],
        estcr[1], estcr[2], estcr[3], estcr[4], estcr[5],
        eprice_T, beta_T[1], beta_T[2]
    )))
end

df = DataFrame(rows)
CSV.write("codefinal.csv", df)


