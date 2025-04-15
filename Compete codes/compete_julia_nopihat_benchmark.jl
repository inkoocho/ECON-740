using Random, Distributions, Plots, LaTeXStrings

alpha=0.96;  #universal benchmark
rho=0.99;    #universal benchmark

sigmab=0.00005;   #pihat experiment
sigmap=0.2;      #pihat experiment
sigmaf=0.0000000001;    #pihat experiment

pihatgain=0.00000001;
delta=0.8;
gam = (1-alpha);
sce = delta/(1-alpha*rho);
Nsim=100;
Titer=10000; 

Vbar=sqrt(sigmap*sigmab*(1-rho^2)/sigmaf);
gain=sqrt(sigmab*sigmaf/(sigmap*(1-rho^2)));

result = zeros(Nsim,2);
price = zeros(Nsim,Titer);
beta0 = zeros(Nsim,Titer);
beta1 = zeros(Nsim,Titer);
pi1 = zeros(Nsim,Titer);
pihat1 = zeros(Nsim,Titer);

for j = 1:Nsim
    V0=zeros(Titer); V1=zeros(Titer); 
    beta0v=zeros(2,Titer); 
    f=zeros(Titer); A1=zeros(Titer); A0=zeros(Titer);
    
    V0[1] = 1.0*Vbar;
    V1[1] = Vbar;

    price[j,1] = gam/(1-alpha)
    beta0[j,1] = gam/(1-alpha) + .01*randn()*sqrt(V0[1])
    beta1[j,1] = price[j,1] + .01*randn()*sqrt(V1[1])
    pi1[j,1] = 0.5;    # universal benchmark 
    #pihat1[j,1] = 0.5;   # comment this line out to assume pihat=0
    #pi1[j,1] = 0; pihat1[1]=0;   #pihat experiment


    vp=sqrt(sigmap)*randn(Titer)
    vf=sqrt(sigmaf)*randn(Titer)
    g=sqrt(sigmab*sigmaf)/sqrt(sigmap*(1-rho^2))
    f[1] = vf[1];
    A1[1]=1; A0[1]=1;
    beta0v[:,1] = [gam/(1-alpha); beta0[j,1]];

    r=zeros(Titer);

    for i=2:Titer             

        f[i] = rho*f[i-1] + vf[i]

        res1 = price[j,i-1] - beta1[j,i-1]
        res0 = price[j,i-1] - (pihat1[j,i-1]*beta1[j,i-1] + (1-pihat1[j,i-1])*(beta0[j,i-1] + sce*f[i-1]))
        #res0 = price[j,i-1] - (pi1[j,i-1]*beta1[j,i-1] + (1-pi1[j,i-1])*(beta0[j,i-1] + sce*f[i-1]))

        beta1[j,i] = beta1[j,i-1] + (V1[i-1]/(sigmap + V1[i-1]))*res1
        
        V1[i] = V1[i-1] - ((V1[i-1])^2)/(sigmap + V1[i-1]) + sigmab
        V0[i] = V0[i-1] - ((V0[i-1])^2)/(sigmap + V0[i-1])

        A1[i] = (exp(-.5*(res1^2)/(sigmap + V1[i-1])))/sqrt(2*3.14159*(sigmap + V1[i-1]))

        M0gaincp = V0[i-1]/(sigmap + V0[i-1])
        
        beta0[j,i] = beta0[j,i-1] + M0gaincp*res0

        A0[i] = (exp(-.5*(res0^2)/(sigmap + V0[i-1]))) / sqrt(2*3.14159*(sigmap + V0[i-1]))

        r[i] = r[i-1] + log(A1[i]/A0[i])
        pi1[j,i] = 1/(1 + exp(-r[i]))
        #pihat1[j,i] = pihat1[j,i-1] + (1/(50+i))*(pi1[j,i] - pihat1[j,i-1])     # comment this line to assume pihat=0
        #pihat1[j,i] = pihat1[j,i-1] + (1/i)*(pi1[j,i] - pihat1[j,i-1])     # comment this line to assume pihat=0
        #pihat1[j,i] = pihat1[j,i-1] + pihatgain*(pi1[j,i] - pihat1[j,i-1])     # comment this line to assume pihat=0

        price[j,i] = gam + delta*f[i] + alpha*((pi1[j,i] + (1-pi1[j,i])*pihat1[j,i])*beta1[j,i] + (1-pi1[j,i])*(1-pihat1[j,i])*beta0[j,i] +
                    (1-pi1[j,i])*(1-pihat1[j,i])*rho*sce*f[i]) + vp[i]
        #price[j,i] = gam + delta*f[i] + alpha*((pi1[j,i] + (1-pi1[j,i])*pi1[j,i])*beta1[j,i] + (1-pi1[j,i])*(1-pi1[j,i])*beta0[j,i] +
        #            (1-pi1[j,i])*(1-pi1[j,i])*rho*sce*f[i]) + vp[i]             
    end
    std_price = sqrt(var(price[j,:]));
    result[j,:]=[pi1[j,Titer], std_price];
end

av = mean.(eachcol(result))
vv = var.(eachcol(result))
println("Average of pi: ", av[1] ,  "     variance of pi: ", vv[1])

#pi_ubar = (4*alpha*rho - 1 + sqrt(1 - 16*alpha*rho*(1-alpha*rho)))/(4*alpha*rho)
#println("pi_ubar: ", pi_ubar)

pp = plot(layout=4)
plot!(price[Nsim,:], title="Price", legend=false, subplot=1)
plot!(pi1[Nsim,:], title="\$ \\pi \$", legend=false,subplot=2)
plot!(beta0[Nsim,:], title="\$ \\beta_0 \$", legend=false,subplot=3)
plot!(beta1[Nsim,:], title="\$ \\beta_1 \$", legend=false,subplot=4)
display(pp)
