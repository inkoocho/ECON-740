##
##    Written by Jung Jae Kim on December 29, 2023. Translated from gresham4_matlab_function.m file
## Modified by Cristhian Rosales-Castillo to estimate the Gresham model using Cobweb model 
## Last modified: November 24, 2025
## Class: Learning and AI in Economics
## Professor: In-Koo Cho

function gresham_julia_function_ehsw(Titer,alpha,rho,sigmap,sigmaf,sigmab)

    #N0 = 0;
    #Nsim = 4000;
    #Titer = 4000;
    #alpha = 0.96;
    #rho = 0.95;
    #sigmap = 1.0;
    #sigmaf = 1.0;
    #sigmab = 0.000015;

    delta = (1 - alpha); # Changed
    sce = delta / (1 - alpha); # Changed
    
    V0 = zeros(Titer);
    V1 = zeros(Titer);
    V0[2] = 0.2; # changed
    V1[2] = 0.2; # changed
    
    beta0 = zeros(Titer);
    beta1 = zeros(Titer);
    pi = zeros(Titer);
    r= zeros(Titer);
    beta0[2] = sce + randn() * sqrt(V0[1]); # changed
    beta1[2] = sce + randn() * sqrt(V1[1]); # changed
    pi[2] = 0.1; # changed
    r[2] = 0.0; # changed
    
    vp = sqrt(sigmap) * randn(Titer);
    vf = sqrt(sigmaf) * randn(Titer);
    #g = sqrt(sigmab * sigmaf) / sqrt(sigmap * (1 - rho^9));
    #varbeta = (g * sigmap) / (sigmaf * (1 - alpha) * (2 - g * (1 - alpha)));
    #varf = sigmaf / (1 - rho^2);
    #varratio = 100 * (varbeta * varf) / (varf * sce^2 + sigmap);
    
    price = zeros(Titer);
    f = zeros(Titer);
    A1 = zeros(Titer);
    A0 = zeros(Titer);
    price[2] = 0.0; # changed
    f[1] = 0.0;
    A1[1] = 1.0;
    A0[1] = 1.0;
    
    MSE1 = zeros(Titer-1);
    MSE0 = zeros(Titer-1);
    
        for i = 3:Titer
            f[i-1] = rho * f[i-2] + vf[i-1] # Changed
    
            res1 = price[i-1] - beta1[i-1] * f[i-2] # changed
            res0 = price[i-1] - beta0[i-1] * f[i-2] # changed
            MSE1[i-1] = sigmap + V1[i-1] * f[i-2]^2 + f[i-2]^2 * V1[i-1] * (1 - alpha * pi[i-1])^2 # changed
            MSE0[i-1] = sigmap + (alpha * pi[i-1])^2 * V1[i-1] * f[i-2]^2 + f[i-2]^2 * V0[i-1] *       
                        (1 - alpha * (1 - pi[i-1]))^2 + V0[i-1] * f[i-2]^2                                 # changed

    
            beta1[i] = beta1[i-1] + (V1[i-1] / MSE1[i-1]) * f[i-2] * res1 #changed
            beta0[i] = beta0[i-1] + (V0[i-1] / MSE0[i-1]) * f[i-2] * res0 #changed
            V1[i] = V1[i-1] - (f[i-2] * V1[i-1])^2 / MSE1[i-1] + sigmab # changed
            V0[i] = V0[i-1] - (f[i-2] * V0[i-1])^2 / MSE0[i-1] + 0.0 # changed
    
            A1[i] = exp(-0.5 * res1^2 / MSE1[i-1]) / sqrt(2 * 3.14159 * MSE1[i-1])
            A0[i] = exp(-0.5 * res0^2 / MSE0[i-1]) / sqrt(2 * 3.14159 * MSE0[i-1])
            r[i] = r[i-1] + log(A1[i] / A0[i])
            pi[i] = 1 / (1 + exp(-r[i]))
            price[i] = (delta + alpha * (pi[i] * beta1[i] + (1 - pi[i]) * beta0[i])) * f[i-1] + vp[i] - alpha * (pi[i] * beta1[i] + (1 - pi[i]) * beta0[i]) * vf[i-1]  #changed
        end
        
    gresham_path_ehsw = hcat(price, pi, beta0, beta1);    
        
    return gresham_path_ehsw    
    end