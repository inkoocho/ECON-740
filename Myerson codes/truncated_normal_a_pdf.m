function pdf = truncated_normal_a_pdf ( x, mu, sigma, a )
    
    if ( x < a )
  
        pdf = 0.0;
    
    else
  
        alpha = ( a - mu ) / sigma;
        xi = ( x - mu ) / sigma;

        alpha_cdf = normcdf ( alpha );
        xi_pdf = normpdf ( xi );

        pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma;

    end
end
