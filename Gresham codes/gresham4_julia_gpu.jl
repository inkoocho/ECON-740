##
##    Written by Jung Jae Kim on March 24, 2025. GPU execution code of gresham4_julia_function.jl
##

using CUDA, Random, Distributions, Plots

function gresham_kernel!(price, pi, beta0, beta1, vf_vals, vp_vals, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute thread index
    if i > Nsim return end  # Ensure within range

    # Initialize per-simulation variables
    V0 = 0.2
    V1 = 0.2
    beta0_i = sce + vp_vals[i, 1] * sqrt(V0)
    beta1_i = sce + vp_vals[i, 2] * sqrt(V1)
    pi_i = 0.1
    r = 0.0
    price_i = 0.0
    f = 0.0

    g = sqrt(sigmab * sigmaf) / sqrt(sigmap * (1 - rho^9))
    varbeta = (g * sigmap) / (sigmaf * (1 - alpha) * (2 - g * (1 - alpha)))

    for t = 2:Titer
        vf = vf_vals[i, t]
        vp = vp_vals[i, t]
        
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
        pi[i, t] = pi_i
        beta0[i, t] = beta0_i
        beta1[i, t] = beta1_i
    end
end

# **Setup Parameters**
Nsim = 4000  # Number of simulations
Titer = 6000  # Number of time steps per simulation
alpha = 0.96
rho = 0.95
sigmap = 1.0
sigmaf = 1.0
sigmab = 0.0005
delta = (1 - rho * alpha)
sce = delta / (1 - alpha * rho)

# **Allocate GPU Memory**
price_gpu = CuArray(zeros(Nsim, Titer))
pi_gpu = CuArray(zeros(Nsim, Titer))
beta0_gpu = CuArray(zeros(Nsim, Titer))
beta1_gpu = CuArray(zeros(Nsim, Titer))

# **Generate Random Numbers on GPU Before Kernel Execution**
vf_gpu = CuArray(randn(Nsim, Titer))  # Pre-generate vf values on GPU
vp_gpu = CuArray(randn(Nsim, Titer));  # Pre-generate vp values on GPU

# **Launch Kernel**
threads_per_block = 256
blocks = cld(Nsim, threads_per_block)
@time @cuda threads=threads_per_block blocks=blocks gresham_kernel!(price_gpu, pi_gpu, beta0_gpu, beta1_gpu, vf_gpu, vp_gpu, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab, delta, sce)

# **Move Results Back to CPU**
price = Array(price_gpu)
pi = Array(pi_gpu)
beta0 = Array(beta0_gpu)
beta1 = Array(beta1_gpu);

# **Plot Results**
plot1 = plot(price[1, :], label="Price")
plot2 = plot(pi[1, :], ylim=(0, 1), label="Pi")
plot3 = plot(beta0[1, :], label="Beta0")
plot4 = plot(beta1[1, :], label="Beta1")
plot(plot1, plot2, plot3, plot4, layout=(2, 2))