##
##    Written by Jung Jae Kim on March 24, 2025. GPU execution code of gresham4_julia_function.jl
##

# Maximum bound of the simulation
# Nsim = 40000, Titer = 10000

using CUDA, Random, Distributions

# SLURM_PROCID gives the task index: 0, 1, 2, 3
# Get task ID (0, 1, 2, 3) from SLURM
task_id = parse(Int, get(ENV, "SLURM_PROCID", "0"))
dev = first(CUDA.devices())
gpu_uuid = CUDA.uuid(dev)
# No need to set CUDA.device! — only one GPU is visible
println("Running SLURM task $task_id using GPU $gpu_uuid")

# **Setup Parameters**
Nsim = 10000  # Number of simulations
Titer = 10000  # Number of time steps per simulation
alpha = 0.96
rho = 0.95
sigmap = 1.0
sigmaf = 1.0
sigmab = 0.0005
delta = (1 - rho * alpha)
sce = delta / (1 - alpha * rho)


function gresham_kernel!(price, pi1, beta0, beta1, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute thread index
    if i > Nsim return end  # Ensure within range

    # Initialize per-simulation variables
    V0 = 0.2
    V1 = 0.2
    beta0_i = sce + randn() * sqrt(V0)
    beta1_i = sce + randn() * sqrt(V1)
    pi_i = 0.1
    r = 0.0
    price_i = 0.0
    f = 0.0

    g = sqrt(sigmab * sigmaf) / sqrt(sigmap * (1 - rho^9))
    varbeta = (g * sigmap) / (sigmaf * (1 - alpha) * (2 - g * (1 - alpha)))

    for t = 2:Titer
        vf = sqrt(sigmap) * randn()
        vp = sqrt(sigmap) * randn()
        
        f = rho * f + vf
        res1 = price_i - beta1_i * f
        res0 = price_i - beta0_i * f
        
        MSE1 = sigmap + V1 * f^2 + f^2 * V1 * (1 - alpha * rho * pi_i)^2
        MSE0 = sigmap + (alpha * rho * pi_i)^2 * V1 * f^2 + f^2 * V0 * (1 - alpha * rho * (1 - pi_i))^2 + V0 * f^2

        beta1_i += (V1 / MSE1) * f * res1
        beta0_i += (V0 / MSE0) * f * res0
        V1 = V1 - (f * V1)^2 / MSE1 + sigmab
        V0 = V0 - (f * V0)^2 / MSE0

        A1 = exp(-0.5 * res1^2 / MSE1) / sqrt(2 * π * MSE1)
        A0 = exp(-0.5 * res0^2 / MSE0) / sqrt(2 * π * MSE0)
        r += log(A1 / A0)
        pi_i = 1 / (1 + exp(-r))
        price_i = (delta + alpha * rho * (pi_i * beta1_i + (1 - pi_i) * beta0_i)) * f + vp
        
        # Store results in global arrays
        price[i, t] = price_i
        pi1[i, t] = pi_i
        beta0[i, t] = beta0_i
        beta1[i, t] = beta1_i
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
@time @cuda threads=threads_per_block blocks=blocks gresham_kernel!(price_gpu, pi_gpu, beta0_gpu, beta1_gpu, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce)

# **Move Results Back to CPU**
price = Array(price_gpu)
pi1 = Array(pi_gpu)
beta0 = Array(beta0_gpu)
beta1 = Array(beta1_gpu);

using CSV, DataFrames
# **Save Results to CSV**
df = DataFrame(price=price[:,Titer], pi1=pi1[:,Titer], beta0=beta0[:,Titer], beta1=beta1[:,Titer])
CSV.write("gresham_results$(task_id).csv", df)
