using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere begin
    alpha=0.96;
    rho=0.99;
    #sigmab=0.05;
    sigmabgrid=0.01;
    sigmabmin=0.000001; sigmabmax= 0.01;
    sigmabinc=sigmabgrid*(sigmabmax-sigmabmin);
    sigmab=sigmabmin:sigmabinc:sigmabmax;
    Numb=length(sigmab);

    sigmap=0.2;

    sigmafgrid=0.01;
    sigmafmin =0.0000000001; sigmafmax=0.000001;
    sigmafinc=sigmafgrid*(sigmafmax-sigmafmin);
    sigmaf=sigmafmin: sigmafinc: sigmafmax;
    Numf=length(sigmaf);

    delta = 0.8;
    gam = (1-alpha);
    sce = delta/(1-alpha*rho);
    Nsim=400;
    Titer=10000;
    #A=Nump+1; AA=A*Numf+Nump+1;
    A=Numb+1; AA=A*Numf+Numb+1; 
    #xt=zeros(Nsim,Numf,Nump);
    xt=zeros(Nsim,Numf,Numb);
    #for j in 1:Nsim, k in 1:Numf, m in 1:Nump
    #    xt[j,k,m]=j*AA+k*A +m
    #end
    for j in 1:Nsim, k in 1:Numf, m in 1:Numb
        xt[j,k,m]=j*AA+k*A +m
    end
end

results = @time @sync @distributed (vcat) for lll in xt
    j=Int(div(lll,AA))
    ll=Int(mod(lll,AA))
    k =Int(div(ll,A))
    m =Int(mod(ll,A))

    V0=zeros(Titer); V1=zeros(Titer); price=zeros(Titer); beta0=zeros(Titer); beta0cp=zeros(Titer); beta0tv=zeros(Titer); beta1=zeros(Titer);
    Vbar=zeros(Titer); gain=zeros(Titer); pi1=zeros(Titer); pihat1=zeros(Titer); beta0v=zeros(2,Titer); 
    f=zeros(Titer); vf=zeros(Titer); A1=zeros(Titer); A0=zeros(Titer);     

    #Vbar=sqrt(sigmap[m]*sigmab*(1-rho^2)/sigmaf[k]);
    Vbar=sqrt(sigmap*sigmab[m]*(1-rho^2)/sigmaf[k]);
    #gain=sqrt(sigmab*sigmaf[k]/(sigmap[m]*(1-rho^2)));
    gain=sqrt(sigmab[m]*sigmaf[k]/(sigmap*(1-rho^2)));

    V0[1] = 1.0*Vbar;
    V1[1] = Vbar;

    price[1] = gam/(1-alpha)
    beta0[1] = gam/(1-alpha) + .01*randn()*sqrt(V0[1])
    #beta0cp[1] = gam/(1-alpha) + .01*randn()*sqrt(V0[1])
    #beta0tv[1] = beta0cp[1]
    beta1[1] = price[1] + .01*randn()*sqrt(V1[1])
    pi1[1] = 0.5; #r = 0.0;
    #pihat1[1] = 0.5;

    #vp=sqrt(sigmap[m])*randn(Titer)
    vp=sqrt(sigmap)*randn(Titer)
    vf=sqrt(sigmaf[k])*randn(Titer)
    #g=sqrt(sigmab*sigmaf[k])/sqrt(sigmap[m]*(1-rho^2))
    g=sqrt(sigmab[m]*sigmaf[k])/sqrt(sigmap*(1-rho^2))
    f[1] = vf[1];
    A1[1]=1; A0[1]=1;
    #phi[:,1] = [1; f[1]];
    #beta0v[:,1] = [gam/(1-alpha); beta0cp[1]];
    beta0v[:,1] = [gam/(1-alpha); beta0[1]];

    r=zeros(Titer);

    for i=2:Titer             

        f[i] = rho*f[i-1] + vf[i]

        res1 = price[i-1] - beta1[i-1]
        res0 = price[i-1] - (pihat1[i-1]*beta1[i-1] + (1-pihat1[i-1])*(beta0[i-1] + sce*f[i-1]))
        #res0 = price[i-1] - (pihat[i-1]*beta1[i-1] + (1-pi1[i-1])*(beta0[i-1] + sce*f[i-1]))

        beta1[i] = beta1[i-1] + (V1[i-1]/(sigmap + V1[i-1]))*res1
        
        V1[i] = V1[i-1] - ((V1[i-1])^2)/(sigmap + V1[i-1]) + sigmab[m]
        V0[i] = V0[i-1] - ((V0[i-1])^2)/(sigmap + V0[i-1])

        A1[i] = (exp(-.5*(res1^2)/(sigmap + V1[i-1])))/sqrt(2*3.14159*(sigmap + V1[i-1]))

        M0gaincp = V0[i-1]/(sigmap + V0[i-1])
        
        beta0[i] = beta0[i-1] + M0gaincp*res0

        A0[i] = (exp(-.5*(res0^2)/(sigmap + V0[i-1]))) / sqrt(2*3.14159*(sigmap + V0[i-1]))

        r[i] = r[i-1] + log(A1[i]/A0[i])
        pi1[i] = 1/(1 + exp(-r[i]))
        #pihat1[i] = pihat1[i-1] + (1/i)*(pi1[i] - pihat1[i-1])

        price[i] = gam + delta*f[i] + alpha*((pi1[i] + (1-pi1[i])*pihat1[i])*beta1[i] + (1-pi1[i])*(1-pihat1[i])*beta0[i] +
                    (1-pi1[i])*(1-pihat1[i])*rho*sce*f[i]) + vp[i]
    end

    compete_path = hcat(j, k, m, pi1[Titer]);
end

#using JLD2
#@save "compete_result.jld"

using DataFrames, CSV, Dates

#result_arrange = zeros(Numf,Nump,Nsim)
#for i=1:Numf*Nump*Nsim
#    j=Integer(results[i,:][1])
#    k=Integer(results[i,:][2])
#    m=Integer(results[i,:][3])
#    result_arrange[k,m,j]=results[i,:][4]
#end

result_arrange = zeros(Numf,Numb,Nsim)
for i=1:Numf*Numb*Nsim
    j=Integer(results[i,:][1])
    k=Integer(results[i,:][2])
    m=Integer(results[i,:][3])
    result_arrange[k,m,j]=results[i,:][4]
end

#av=zeros(Numf,Nump);
#for k=1:Numf, m=1:Nump
#   av[k,m]=sum(result_arrange[k,m,:])/length(result_arrange[k,m,:]);
#end
av=zeros(Numf,Numb);
for k=1:Numf, m=1:Numb
   av[k,m]=sum(result_arrange[k,m,:])/length(result_arrange[k,m,:]);
end

#df = DataFrame(avg_Pi = vec(av),sigmaf = repeat(sigmaf,outer=Nump), sigmap = repeat(sigmap,inner=Numf), sigmab=sigmab )
df = DataFrame(avg_Pi = vec(av),sigmaf = repeat(sigmaf,outer=Numb), sigmab = repeat(sigmab,inner=Numf), sigmap=sigmap, alpha=alpha )
nameofoutput="result"; nameofformat=".csv";nameoffile=nameofoutput*Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")*nameofformat
CSV.write(nameoffile,df)