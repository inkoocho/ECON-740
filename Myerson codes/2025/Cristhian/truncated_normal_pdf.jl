function truncated_normal_pdf(x, mu, sigma, a, b)
    d = truncated(Normal(mu, sigma), a, b)
    return pdf(d, x)
end