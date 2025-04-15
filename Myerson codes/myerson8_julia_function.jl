function myerson8_julia_function(mu, sigma, Nbuyer, a, epsilon, T, CRsample)

    # Helper functions
    truncated_normal_a_pdf(x, mu, sigma, a) = x < a ? 0.0 :
        begin
            alpha = (a - mu) / sigma
            xi = (x - mu) / sigma
            alpha_cdf = cdf(Normal(0,1), alpha)
            xi_pdf = pdf(Normal(0,1), xi)
            xi_pdf / ((1.0 - alpha_cdf) * sigma)
        end

    truncated_normal_a_cdf(x, mu, sigma, a) = x < a ? 0.0 :
        begin
            alpha = (a - mu) / sigma
            xi = (x - mu) / sigma
            alpha_cdf = cdf(Normal(0,1), alpha)
            xi_cdf = cdf(Normal(0,1), xi)
            (xi_cdf - alpha_cdf) / (1.0 - alpha_cdf)
        end

    truncated_normal_a_sample(mu, sigma, a, n) = begin
        alpha = (a - mu) / sigma
        alpha_cdf = cdf(Normal(0,1), alpha)
        u = rand(n)
        xi_cdf = alpha_cdf .+ u .* (1.0 - alpha_cdf)
        xi = quantile.(Normal(0,1), xi_cdf)
        mu .+ sigma .* xi
    end

    hazard(x, mu, sigma, a) = (1.0 - truncated_normal_a_cdf(x, mu, sigma, a)) / truncated_normal_a_pdf(x, mu, sigma, a)

    hazardrate(x) = hazard(x, mu, sigma, mu) - x
    pstar = find_zero(hazardrate, (mu, 3*mu), Bisection(); verbose=false)
    fval = hazardrate(pstar)

    if abs(fval) > 1e-3
        @warn "Optimal pricing is not found for j = $j, fval = $fval"
    end

    pistar = pstar * (1 - truncated_normal_a_cdf(pstar, mu, sigma, mu))

    beta = zeros(2, T+1)
    beta[:, 1] .= [15.0, -0.5]

    R = zeros(2, 2, T)
    R[:, :, 1] = Matrix{Float64}(I, 2, 2)

    e = 2 * epsilon .* (rand(Binomial(1, 0.5), T) .- 0.5)
    eprice = zeros(T)
    p = zeros(T)
    Fp = zeros(T)
    Dagg = zeros(T)
    q = zeros(T)
    pi = zeros(T)

    for iter in 1:T
        eprice[iter] = -beta[1, iter] / (2 * beta[2, iter])
        p[iter] = eprice[iter] + e[iter]
        Fp[iter] = truncated_normal_a_cdf(p[iter], mu, sigma, mu)
        Dagg[iter] = rand(Binomial(Nbuyer, 1 - Fp[iter]))
        q[iter] = Dagg[iter] / Nbuyer
        pi[iter] = q[iter] * p[iter]

        R[:, :, iter] = [1 eprice[iter]; eprice[iter] eprice[iter]^2 + epsilon^2]

        if q[iter] == 0
            beta[:, iter+1] = beta[:, iter] + a * ([-epsilon, 0.0] .- beta[:, iter])
        elseif q[iter] == 1
            beta[:, iter+1] = beta[:, iter] + a * ([epsilon, 0.0] .- beta[:, iter])
        else
            beta[:, iter+1] = beta[:, iter] + a * (R[:, :, iter] \ ([1.0, p[iter]] * (q[iter] - dot([1.0, p[iter]], beta[:, iter]))))
        end
    end

    CLferror = eprice[T] - pstar
    forecastp = eprice
    realizedp = p
    realizedq = q
    estimated = beta

    CRferror = zeros(CRsample)
    estimatedcr = zeros(CRsample)
    for isample in 1:CRsample
        Obvalue = 2 * T * isample
        values = truncated_normal_a_sample(mu, sigma, mu, Obvalue)
        qvec = (1:Obvalue) ./ Obvalue
        pcr = sort(values, rev=true)
        rev = qvec .* pcr
        M, I = findmax(rev)
        CRferror[isample] = pstar - pcr[I]
        estimatedcr[isample] = pcr[I]
    end

    hcat(CLferror, CRferror')
end