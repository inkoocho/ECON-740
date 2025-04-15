function cdf = truncated_normal_a_cdf ( x, mu, sigma, a )
    if ( x < a )
  
        cdf = 0.0;
    
    else
  
        alpha = ( a - mu ) / sigma;
        xi = ( x - mu ) / sigma;

        alpha_cdf = normcdf ( alpha );
        xi_cdf = normcdf ( xi );

        cdf = ( xi_cdf - alpha_cdf ) / ( 1.0 - alpha_cdf );
    
    end
end
