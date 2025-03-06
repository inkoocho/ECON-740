##
##    Written by Jung Jae Kim on December 29, 2023. Translated from gresham4_matlab_forparallel.m file
##

function gresham4_julia_function(Titer,alpha,rho,sigmap,sigmaf,sigmab)

    #N0 = 0;
    #Nsim = 4000;
    #Titer = 4000;
    #alpha = 0.96;
    #rho = 0.95;
    #sigmap = 1.0;
    #sigmaf = 1.0;
    #sigmab = 0.000015;

    delta = (1 - rho * alpha);
    sce = delta / (1 - alpha * rho);
    
    V0 = zeros(Titer);
    V1 = zeros(Titer);
    V0[1] = 0.2;
    V1[1] = 0.2;
    
    beta0 = zeros(Titer);
    beta1 = zeros(Titer);
    pi = zeros(Titer);
    r= zeros(Titer);
    beta0[1] = sce + randn() * sqrt(V0[1]);
    beta1[1] = sce + randn() * sqrt(V1[1]);
    pi[1] = 0.1;
    r[1] = 0.0;
    
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
    price[1] = 0.0;
    f[1] = 0.0;
    A1[1] = 1.0;
    A0[1] = 1.0;
    
    MSE1 = zeros(Titer-1);
    MSE0 = zeros(Titer-1);
    
        for i = 2:Titer
            f[i] = rho * f[i-1] + vf[i]
    
            res1 = price[i-1] - beta1[i-1] * f[i-1]
            res0 = price[i-1] - beta0[i-1] * f[i-1]
            MSE1[i-1] = sigmap + V1[i-1] * f[i-1]^2 + f[i-1] * f[i-1] * V1[i-1] * (1 - alpha * rho * pi[i-1])^2
            MSE0[i-1] = sigmap + (alpha * rho * pi[i-1])^2 * V1[i-1] * f[i-1]^2 + f[i-1]^2 * V0[i-1] *
                        (1 - alpha * rho * (1 - pi[i-1]))^2 + V0[i-1] * f[i-1]^2
    
            beta1[i] = beta1[i-1] + (V1[i-1] / MSE1[i-1]) * f[i-1] * res1
            beta0[i] = beta0[i-1] + (V0[i-1] / MSE0[i-1]) * f[i-1] * res0
            V1[i] = V1[i-1] - (f[i-1] * V1[i-1])^2 / MSE1[i-1] + sigmab
            V0[i] = V0[i-1] - (f[i-1] * V0[i-1])^2 / MSE0[i-1] + 0.0
    
            A1[i] = exp(-0.5 * res1^2 / MSE1[i-1]) / sqrt(2 * 3.14159 * MSE1[i-1])
            A0[i] = exp(-0.5 * res0^2 / MSE0[i-1]) / sqrt(2 * 3.14159 * MSE0[i-1])
            r[i] = r[i-1] + log(A1[i] / A0[i])
            pi[i] = 1 / (1 + exp(-r[i]))
            price[i] = (delta + alpha * rho * (pi[i] * beta1[i] + (1 - pi[i]) * beta0[i])) * f[i] + vp[i]  
        end
        
    gresham_path = hcat(price, pi, beta0, beta1);    

    return gresham_path    
    end