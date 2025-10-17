function myerson8_function(mu, sigma, N, a, epsilon, T, CRsample)
    
    hazardrate(x) = hazard_rate(truncated(Normal(mu, sigma), mu, Inf), x) - x

    pstar = find_zero(hazardrate, (mu, 3*mu), Bisection())
    fval = hazard_rate(truncated(Normal(mu, sigma), mu, Inf), pstar)

    if abs(fval) > 1e-3
        @warn "Optimal pricing is not found given current parametrization."
    end

    pistar = pstar*(1 - truncated_normal_cdf(pstar, mu, sigma, mu, Inf))

    beta0 = zeros(T+1)
    beta1 = zeros(T+1)
    beta0[1] = 15
    beta1[1] = -0.5
    R = [zeros(2,2) for _ in 1:T]
    R[1] = [1.0 0.0; 0.0 1.0]

    e = 2 * epsilon * (rand(Binomial(1, 0.5), T) .- 0.5)
    eprice  = zeros(T)
    p = zeros(T)
    Fp = zeros(T)
    Dagg = zeros(T)
    q = zeros(T)
    pi = zeros(T)

    for iter in 1:T

        eprice[iter] = -beta0[iter] / (2 * beta1[iter])
        p[iter] = eprice[iter] + e[iter]
        Fp[iter] = truncated_normal_cdf(p[iter], mu, sigma, mu, Inf)
        Dagg[iter] = rand(Binomial(N, 1 - Fp[iter]))
        q[iter] = Dagg[iter] / N
        pi[iter] = p[iter] * q[iter]

        R[iter] = [1.0 eprice[iter]; eprice[iter] eprice[iter]^2+epsilon^2]

        if q[iter] == 0
            update_vect = [beta0[iter]; beta1[iter]] + a * ([-epsilon; 0] - [beta0[iter]; beta1[iter]])
            beta0[iter+1] = update_vect[1]
            beta1[iter+1] = update_vect[2]
        elseif q[iter] == 1
            update_vect = [beta0[iter]; beta1[iter]] + a * ([epsilon; 0] - [beta0[iter]; beta1[iter]])
            beta0[iter+1] = update_vect[1]
            beta1[iter+1] = update_vect[2]
        else
            q_error = q[iter] -(beta0[iter] + beta1[iter] * p[iter])
            update_vect = [beta0[iter]; beta1[iter]] + a * (R[iter] \ ([1 ; p[iter]] * q_error))
            beta0[iter+1] = update_vect[1]
            beta1[iter+1] = update_vect[2]
        end
    end

    CLferror = eprice[T] - pstar
    forecastp = eprice
    realizedp = p
    realizedq = q
    estimated_beta0 = beta0
    estimated_beta1 = beta1

    CRferror = zeros(CRsample)
    estimatedcr = zeros(CRsample)

    for isample in 1:CRsample
        Obvalue = 2 * T * isample
        value = truncated_normal_sample(mu, sigma, mu, Inf, Obvalue)
        q = (1:Obvalue) ./ Obvalue
        pcr = sort(value, rev=true)
        rev = q .* pcr
        M, I = findmax(rev)
        CRferror[isample] = pstar - pcr[I]
        estimatedcr[isample] = pcr[I]
    end
    return CLferror, CRferror
end
