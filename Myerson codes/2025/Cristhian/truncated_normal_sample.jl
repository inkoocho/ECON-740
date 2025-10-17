function truncated_normal_sample(mu, sigma, a, b, n)
    d = truncated(Normal(mu, sigma), a, b)
    return rand(d, n)
end