function hazard_rate(dist, x)
    F = cdf(dist, x)
    f = pdf(dist, x)
    return (1 - F) / f
end