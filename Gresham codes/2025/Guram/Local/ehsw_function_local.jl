##
##    Modified by Guram Lobzhanidze on November 17, 2025. Original author: Jung Jae Kim
##

function ehsw_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)

    #N0 = 0;
    #Nsim = 4000;
    #Titer = 4000;
    #alpha = 0.96;
    #rho = 0.95;
    #sigmap = 1.0;
    #sigmaf = 1.0;
    #sigmab = 0.000015;

    delta = (1 - rho * alpha)
    sce   = delta / (1 - alpha * rho)

    V0 = zeros(Titer)
    V1 = zeros(Titer)
    V0[1] = 0.2
    V1[1] = 0.2

    beta0 = zeros(Titer)
    beta1 = zeros(Titer)
    pi    = zeros(Titer)
    r     = zeros(Titer)
    beta0[1] = sce + randn() * sqrt(V0[1])
    beta1[1] = sce + randn() * sqrt(V1[1])
    pi[1]    = 0.1
    r[1]     = 0.0

    vp = sqrt(sigmap) * randn(Titer)
    vf = sqrt(sigmaf) * randn(Titer)

    price = zeros(Titer)
    f     = zeros(Titer)          # fundamental z_t
    zlag  = zeros(Titer)          # NEW: lagged fundamental z_{t-1}
    A1 = zeros(Titer)
    A0 = zeros(Titer)

    price[1] = 0.0
    f[1]     = 0.0
    zlag[1]  = 0.0               # at t=1, z_{0} ≈ 0
    A1[1]    = 1.0
    A0[1]    = 1.0

    MSE1 = zeros(Titer-1)
    MSE0 = zeros(Titer-1)

    for i = 2:Titer
        # AR(1) for the fundamental: z_t
        f[i] = rho * f[i-1] + vf[i]

        # lagged fundamental z_{t-1}
        zlag[i] = f[i-1]

        # --- Bayesian updating based on observation at time i-1 ---
        #   p_{i-1} = β * z_{(i-1)-1} + η_{i-1} = β * zlag[i-1] + η_{i-1}

        res1 = price[i-1] - beta1[i-1] * zlag[i-1]
        res0 = price[i-1] - beta0[i-1] * zlag[i-1]

        MSE1[i-1] = sigmap +
                    V1[i-1] * zlag[i-1]^2 +
                    zlag[i-1]^2 * V1[i-1] * (1 - alpha * rho * pi[i-1])^2

        MSE0[i-1] = sigmap +
                    (alpha * rho * pi[i-1])^2 * V1[i-1] * zlag[i-1]^2 +
                    zlag[i-1]^2 * V0[i-1] * (1 - alpha * rho * (1 - pi[i-1]))^2 +
                    V0[i-1] * zlag[i-1]^2

        beta1[i] = beta1[i-1] + (V1[i-1] / MSE1[i-1]) * zlag[i-1] * res1
        beta0[i] = beta0[i-1] + (V0[i-1] / MSE0[i-1]) * zlag[i-1] * res0

        V1[i] = V1[i-1] - (zlag[i-1] * V1[i-1])^2 / MSE1[i-1] + sigmab
        V0[i] = V0[i-1] - (zlag[i-1] * V0[i-1])^2 / MSE0[i-1] + 0.0

        A1[i] = exp(-0.5 * res1^2 / MSE1[i-1]) / sqrt(2 * 3.14159 * MSE1[i-1])
        A0[i] = exp(-0.5 * res0^2 / MSE0[i-1]) / sqrt(2 * 3.14159 * MSE0[i-1])

        r[i]  = r[i-1] + log(A1[i] / A0[i])
        pi[i] = 1 / (1 + exp(-r[i]))

        # --- NEW equilibrium condition: p_t = β * z_{t-1} + η_t ---
        # here zlag[i] is z_{t-1}
        price[i] = (delta + alpha * rho * (pi[i] * beta1[i] + (1 - pi[i]) * beta0[i])) * zlag[i] + vp[i]
    end

    ehsw_path = hcat(price, pi, beta0, beta1)

    return ehsw_path
end




