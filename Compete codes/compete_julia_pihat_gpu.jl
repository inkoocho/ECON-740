##
##    Written by Jung Jae Kim on March 24, 2025. GPU execution code of compete_julia_pihat_benchmark.jl
##

# Maximum bound of the simulation
# Nsim = 400, Numf = 10, Numb = 10, Titer = 10000

using CUDA, Random, Distributions, Plots

# **Setup Parameters**
Nsim = 400
Titer = 10000

sigmafmin = 0.0000000001
sigmafmax = 0.000001
Numf = 10
sigmaf_cpu = range(sigmafmin, stop=sigmafmax, length=Numf) |> collect
sigmaf_grid = CuArray(sigmaf_cpu)

sigmabmin = 0.000001
sigmabmax = 0.01
Numb = 10 
sigmab_cpu = range(sigmabmin, stop=sigmabmax, length=Numb) |> collect
sigmab_grid = CuArray(sigmab_cpu)

alpha = 0.96
rho = 0.99
sigmap = 0.2

TotalSim = Nsim * Numf * Numb

function compete_kernel!(
    price, pi1, pihat1, beta0, beta1,
    sigmaf_grid, sigmab_grid,
    Nsim, Ns, Nb, Titer,
    alpha, rho, sigmap
)
    id = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    TotalSim = Nsim * Ns * Nb
    if id > TotalSim
        return
    end

    # Unpack simulation indices
    i_sim = (id - 1) % Nsim + 1
    i_sf  = ((id - 1) ÷ Nsim) % Ns + 1
    i_sb  = ((id - 1) ÷ (Nsim * Ns)) % Nb + 1

    # Retrieve sigma values for this thread
    sigmaf = sigmaf_grid[i_sf]
    sigmab = sigmab_grid[i_sb]

    # Initialize parameters
    delta = 0.8
    gam = (1 - alpha)
    sce = delta / (1 - alpha * rho)
    Vbar = sqrt(sigmap * sigmab * (1 - rho^2) / sigmaf)

    V0 = 1.0 * Vbar
    V1 = Vbar

    price_i = gam / (1 - alpha)
    beta0_i = gam / (1 - alpha) + 0.01 * randn() * sqrt(V0)
    beta1_i = price_i + 0.01 * randn() * sqrt(V1)
    pi1_i = 0.5
    pihat1_i = 0.5

    g = sqrt(sigmab * sigmaf) / sqrt(sigmap * (1 - rho^2))
    f = sqrt(sigmaf) * randn()
    A1 = 1.0
    A0 = 1.0
    r = 0.0

    # Save initial values
    price[i_sim, i_sf, i_sb, 1] = price_i
    pi1[i_sim, i_sf, i_sb, 1] = pi1_i
    pihat1[i_sim, i_sf, i_sb, 1] = pihat1_i
    beta0[i_sim, i_sf, i_sb, 1] = beta0_i
    beta1[i_sim, i_sf, i_sb, 1] = beta1_i

    for t = 2:Titer
        vp = sqrt(sigmap) * randn()
        vf = sqrt(sigmaf) * randn()
        
        res1 = price_i - beta1_i
        res0 = price_i - (pihat1_i * beta1_i + (1 - pihat1_i) * (beta0_i + sce * f))
        f = rho * f + vf

        beta1_i += (V1 / (sigmap + V1)) * res1
        A1 = exp(-0.5 * (res1^2) / (sigmap + V1)) / sqrt(2 * 3.14159 * (sigmap + V1))

        gain0 = V0 / (sigmap + V0)
        beta0_i += gain0 * res0
        A0 = exp(-0.5 * (res0^2) / (sigmap + V0)) / sqrt(2 * 3.14159 * (sigmap + V0))

        V1 = V1 - (V1^2) / (sigmap + V1) + sigmab
        V0 = V0 - (V0^2) / (sigmap + V0)

        r += log(A1 / A0)
        pi1_i = 1 / (1 + exp(-r))
        pihat1_i += (1 / (50 + t)) * (pi1_i - pihat1_i)

        price_i = gam + delta * f + alpha * (
            (pi1_i + (1 - pi1_i) * pihat1_i) * beta1_i +
            (1 - pi1_i) * (1 - pihat1_i) * beta0_i +
            (1 - pi1_i) * (1 - pihat1_i) * rho * sce * f
        ) + vp

        # Save results
        price[i_sim, i_sf, i_sb, t] = price_i
        pi1[i_sim, i_sf, i_sb, t] = pi1_i
        pihat1[i_sim, i_sf, i_sb, t] = pihat1_i
        beta0[i_sim, i_sf, i_sb, t] = beta0_i
        beta1[i_sim, i_sf, i_sb, t] = beta1_i
    end
end

# Allocate 4D arrays
price_gpu = CuArray(zeros(Float32, Nsim, Numf, Numb, Titer))
pi1_gpu = CuArray(zeros(Float32, Nsim, Numf, Numb, Titer))
pihat1_gpu = CuArray(zeros(Float32, Nsim, Numf, Numb, Titer))
beta0_gpu = CuArray(zeros(Float32, Nsim, Numf, Numb, Titer))
beta1_gpu = CuArray(zeros(Float32, Nsim, Numf, Numb, Titer))

# Launch kernel
threads = 256
blocks = cld(TotalSim, threads)
@time @cuda threads=threads blocks=blocks compete_kernel!(
    price_gpu, pi1_gpu, pihat1_gpu, beta0_gpu, beta1_gpu,
    sigmaf_grid, sigmab_grid,
    Nsim, Numf, Numb, Titer,
    alpha, rho, sigmap
)

# **Move Results Back to CPU**
price = Array(price_gpu)
pi1 = Array(pi1_gpu)
pihat1 = Array(pihat1_gpu)
beta0 = Array(beta0_gpu)
beta1 = Array(beta1_gpu);

using Plots, LaTeXStrings
# **Plot Results**
pp = plot(layout=4)
plot!(price[1,1,1,:], title="Price", legend=false, subplot=1)
plot!(pihat1[1,1,1,:], title="\$ \\hat{\\pi}\$", legend=false,subplot=2)
plot!(beta0[1,1,1,:], title="\$ \\beta_0 \$", legend=false,subplot=3)
plot!(beta1[1,1,1,:], title="\$ \\beta_1 \$", legend=false,subplot=4)
savefig(pp,"compete_results.png")