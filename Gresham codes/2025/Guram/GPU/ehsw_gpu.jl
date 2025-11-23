##
##    Written by Jung Jae Kim on March 24, 2025. GPU execution code of gresham4_julia_function.jl
##

# Maximum bound of the simulation
# Nsim = 40000, Titer = 10000

using CUDA, Random, Distributions

# **Setup Parameters**
Nsim = 10000  # Number of simulations
Titer = 20000  # Number of time steps per simulation
alpha = 0.96
rho = 0.95
sigmap = 1.0
sigmaf = 1.0
sigmab = 0.0005
delta = (1 - rho * alpha)
sce = delta / (1 - alpha * rho)


function ehsw_kernel!(
    price, pi, beta0, beta1,
    Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce
)
    sim = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if sim > Nsim
        return
    end

    # --- per-simulation scalars (match local init) ---
    V0 = 0.2
    V1 = 0.2
    beta0_i = sce + randn() * sqrt(V0)
    beta1_i = sce + randn() * sqrt(V1)
    pi_i = 0.1
    r_i  = 0.0

    price_prev = 0.0    # p_{t-1}, with p_1 = 0
    f_prev     = 0.0    # z_{t-1}, with z_1 = 0
    zlag_prev  = 0.0    # z_{t-2}, with z_0 ≈ 0

    # store t=1 values (optional but nice for parity)
    price[sim, 1] = price_prev
    pi[sim, 1]    = pi_i
    beta0[sim, 1] = beta0_i
    beta1[sim, 1] = beta1_i

    for t = 2:Titer
        # shocks (match local: vp~N(0,sigmap), vf~N(0,sigmaf))
        vf_t = sqrt(sigmaf) * randn()
        vp_t = sqrt(sigmap) * randn()

        # fundamental AR(1): z_t
        f_cur = rho * f_prev + vf_t

        # lagged fundamental z_{t-1}
        zlag_cur = f_prev

        # residuals use p_{t-1} and z_{t-2} (local uses zlag[t-1])
        res1 = price_prev - beta1_i * zlag_prev
        res0 = price_prev - beta0_i * zlag_prev

        # MSEs use z_{t-2} and beliefs at t-1 (local MSE1[t-1], MSE0[t-1])
        MSE1 = sigmap +
               V1 * zlag_prev^2 +
               zlag_prev^2 * V1 * (1 - alpha * rho * pi_i)^2

        MSE0 = sigmap +
               (alpha * rho * pi_i)^2 * V1 * zlag_prev^2 +
               zlag_prev^2 * V0 * (1 - alpha * rho * (1 - pi_i))^2 +
               V0 * zlag_prev^2

        # Bayesian updating
        beta1_i = beta1_i + (V1 / MSE1) * zlag_prev * res1
        beta0_i = beta0_i + (V0 / MSE0) * zlag_prev * res0

        V1 = V1 - (zlag_prev * V1)^2 / MSE1 + sigmab
        V0 = V0 - (zlag_prev * V0)^2 / MSE0 + 0.0

        # likelihood ratio / posterior
        A1 = exp(-0.5 * res1^2 / MSE1) / sqrt(2 * 3.14159 * MSE1)
        A0 = exp(-0.5 * res0^2 / MSE0) / sqrt(2 * 3.14159 * MSE0)

        r_i  = r_i + log(A1 / A0)
        pi_i = 1 / (1 + exp(-r_i))

        # equilibrium price uses z_{t-1} (local: price[t] depends on zlag[t])
        price_cur = (delta + alpha * rho *
                    (pi_i * beta1_i + (1 - pi_i) * beta0_i)) * zlag_cur + vp_t

        # store
        price[sim, t] = price_cur
        pi[sim, t]    = pi_i
        beta0[sim, t] = beta0_i
        beta1[sim, t] = beta1_i

        # roll forward
        price_prev = price_cur
        zlag_prev  = zlag_cur   # becomes z_{t-1} next step's "z_{t-2}"
        f_prev     = f_cur
    end
end

# **Allocate GPU Memory**
price_gpu = CuArray(zeros(Nsim, Titer))
pi_gpu = CuArray(zeros(Nsim, Titer))
beta0_gpu = CuArray(zeros(Nsim, Titer))
beta1_gpu = CuArray(zeros(Nsim, Titer))

# **Launch Kernel**
threads_per_block = 256
blocks = cld(Nsim, threads_per_block)

@time @cuda threads=threads_per_block blocks=blocks ehsw_kernel!(
    price_gpu, pi_gpu, beta0_gpu, beta1_gpu,
    Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce
)


# **Move Results Back to CPU**
price = Array(price_gpu)
pi = Array(pi_gpu)
beta0 = Array(beta0_gpu)
beta1 = Array(beta1_gpu);



using CSV, DataFrames
# **Save Results to CSV**



df = DataFrame(sim=Int[], time=Int[], price = Float64[], pi = Float64[], beta0 = Float64[], beta1 = Float64[])


for sim_id in 1:Nsim
    for time in 1:Titer
        push!(df, (sim_id, time,
                   price[sim_id, time],
                   pi[sim_id, time],
                   beta0[sim_id, time],
                   beta1[sim_id, time]))
    end
end




CSV.write("ehsw_gpu_results.csv", df)