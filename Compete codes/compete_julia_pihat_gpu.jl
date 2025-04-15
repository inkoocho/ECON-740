using CUDA, Random, Distributions, Plots

function compete_kernel!(price, pi1, pihat1, beta0, beta1, vf_vals, vp_vals, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute thread index
    if i > Nsim return end  # Ensure within range

    # Initialize per-simulation variables
    delta=0.8
    gam = (1-alpha)
    sce = delta/(1-alpha*rho)

    Vbar=sqrt(sigmap*sigmab*(1-rho^2)/sigmaf)

    V0 = 1.0*Vbar
    V1 = Vbar

    price_i = gam/(1-alpha)
    beta0_i = gam/(1-alpha) + .01*vp_vals[i,1]*sqrt(V0)
    beta1_i = price_i + .01*vp_vals[i,2]*sqrt(V1)
    pi1_i = 0.5
    pihat1_i = 0.5

    g=sqrt(sigmab*sigmaf)/sqrt(sigmap*(1-rho^2))
    f = sqrt(sigmaf)*vf_vals[i,1]
    A1=1
    A0=1
    r=0.0

    price[i, 1] = price_i
    pi1[i, 1] = pi1_i
    pihat1[i, 1] = pihat1_i
    beta0[i, 1] = beta0_i
    beta1[i, 1] = beta1_i

    for t = 2:Titer
        vp = sqrt(sigmap)*vf_vals[i, t]
        vf = sqrt(sigmaf)*vp_vals[i, t]
        
        res1 = price_i - beta1_i 
        res0 = price_i - (pihat1_i * beta1_i + (1 - pihat1_i)*(beta0_i + sce*f))
        f = rho * f + vf

        beta1_i = beta1_i + (V1 / (sigmap + V1))*res1

        A1 = (exp(-.5 * (res1^2) / (sigmap + V1))) / sqrt(2*3.14159*(sigmap + V1))

        M0gaincp = V0 / (sigmap + V0)

        beta0_i = beta0_i +M0gaincp*res0

        A0 = (exp(-.5 * (res0^2) / (sigmap + V0))) / sqrt(2*3.14159*(sigmap + V0))

        V1 = V1 - (V1^2) / (sigmap + V1) + sigmab
        V0 = V0 - (V0^2) / (sigmap + V0)

        r = r + log(A1/A0)
        pi1_i = 1 / (1 + exp(-r))
        pihat1_i = pihat1_i + (1/(50+t))*(pi1_i - pihat1_i)

        price_i = gam + delta*f + alpha*((pi1_i + (1 - pi1_i)*pihat1_i)*beta1_i + (1 - pi1_i)*(1 - pihat1_i)*beta0_i +
                (1 - pi1_i)*(1 - pihat1_i)*rho*sce*f) + vp
        
        # Store results in global arrays
        price[i, t] = price_i
        pi1[i, t] = pi1_i
        pihat1[i, t] = pihat1_i
        beta0[i, t] = beta0_i
        beta1[i, t] = beta1_i
    end
end

# **Setup Parameters**
Nsim = 100  # Number of simulations
Titer = 10000  # Number of time steps per simulation
alpha = 0.96
rho = 0.99
sigmap = 0.2
sigmaf = 0.0000000001
sigmab = 0.00005;

# **Allocate GPU Memory**
price_gpu = CuArray(zeros(Nsim, Titer))
pi1_gpu = CuArray(zeros(Nsim, Titer))
pihat1_gpu = CuArray(zeros(Nsim, Titer))
beta0_gpu = CuArray(zeros(Nsim, Titer))
beta1_gpu = CuArray(zeros(Nsim, Titer))

# **Generate Random Numbers on GPU Before Kernel Execution**
vf_gpu = CuArray(randn(Nsim, Titer))  # Pre-generate vf values on GPU
vp_gpu = CuArray(randn(Nsim, Titer));  # Pre-generate vp values on GPU

# **Launch Kernel**
threads_per_block = 256
blocks = cld(Nsim, threads_per_block)
@time @cuda threads=threads_per_block blocks=blocks compete_kernel!(price_gpu, pi1_gpu, pihat1_gpu, beta0_gpu, beta1_gpu, vf_gpu, vp_gpu, Nsim, Titer, alpha, rho, sigmap, sigmaf, sigmab)

# **Move Results Back to CPU**
price = Array(price_gpu)
pi1 = Array(pi1_gpu)
pihat1 = Array(pihat1_gpu)
beta0 = Array(beta0_gpu)
beta1 = Array(beta1_gpu);

using Plots, LaTeXStrings
# **Plot Results**
pp = plot(layout=4)
plot!(price[17,:], title="Price", legend=false, subplot=1)
plot!(pihat1[17,:], title="\$ \\hat{\\pi}\$", legend=false,subplot=2)
plot!(beta0[17,:], title="\$ \\beta_0 \$", legend=false,subplot=3)
plot!(beta1[17,:], title="\$ \\beta_1 \$", legend=false,subplot=4)
display(pp)