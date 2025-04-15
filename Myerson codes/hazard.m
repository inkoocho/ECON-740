function rate = hazard(x, mu, sigma, a)
    
    rate = (1-truncated_normal_a_cdf(x, mu, sigma, a))/truncated_normal_a_pdf(x, mu, sigma, a);
    
end