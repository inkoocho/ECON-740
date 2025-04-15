function x = truncated_normal_a_sample ( mu, sigma, a, n )

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normcdf( alpha );
  beta_cdf = 1.0;

  u = rand ( 1, n );
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
  xi = norminv( xi_cdf );

  x = mu + sigma * xi;

end
