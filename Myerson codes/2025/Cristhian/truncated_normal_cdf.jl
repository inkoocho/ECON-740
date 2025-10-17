function truncated_normal_cdf(x, mu, sigma, a, b)
    d = truncated(Normal(mu, sigma), a, b)
    return cdf(d, x)
end
