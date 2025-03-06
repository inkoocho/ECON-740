function oligopoly_julia_function(N,lambda,sigmawn)

    #N = 20000;
    n_firm =3;

    AA = 1;
    BB = 1;
    CC = 0.7;

    A = fill(AA,n_firm);
    B = fill(BB,n_firm);
    C = fill(CC,n_firm);
    #C = CC .+ 0.05*randn(n_firm) # general C_j

    # NE
    p_nash = AA / (2*BB-CC);
    p_collude=AA /(2*(BB-CC));
    #p_nash = zeros(n_firm);
    #p_collude = zeros(n_firm);
    #for i = 1:n_firm
    #    p_nash[i] = AA / (2*BB - (sum(C)-C[i]))
    #    p_collude[i] = AA / (2*(BB - (sum(C)-C[i])))
    #end

    # set 1he matrix size
    p = zeros(n_firm,N);
    q = zeros(n_firm,N);

    alpha1 = zeros(4,N);
    alpha2 = zeros(4,N);
    alpha3 = zeros(4,N);
    #beta = zeros(n_firm,3,N)

    br = zeros(n_firm,N);
    pavg = zeros(n_firm,N);
    profit = zeros(n_firm,N);
    p12_avg = zeros(1,N);
    p13_avg = zeros(1,N);
    p23_avg = zeros(1,N);

    #sigmawn=0.065;
    wn = sigmawn*randn(n_firm,N);

    R12 = Matrix{Float64}(I, 2, 2);
    R13 = Matrix{Float64}(I, 2, 2);
    R21 = Matrix{Float64}(I, 2, 2);
    R23 = Matrix{Float64}(I, 2, 2);
    R31 = Matrix{Float64}(I, 2, 2);
    R32 = Matrix{Float64}(I, 2, 2);

    # initialize
    alpha1[:,1] = [p_nash 0 p_nash 0]
    alpha2[:,1] = [p_nash 0 p_nash 0]
    alpha3[:,1] = [p_nash 0 p_nash 0]
    #beta[:,:,1] .= [AA; BB; CC]
    #beta[:,:,1] = [A B C]

    #sig = [0.1, 0.1];
    #lambda = 0.0025;

    p[:,1] .= p_nash
    pavg[:,1] = p[:,1];
    p12_avg[1,1] = p_nash;
    p13_avg[1,1] = p_nash;
    p23_avg[1,1] = p_nash;

    for t = 1:N-1
    
        br[1,t] = -(2*(128*AA*BB^2 + 64*BB^2*CC*alpha1[1,t] + 64*BB^2*CC*alpha1[3,t] - 4*AA*CC^2*alpha1[2,t]*alpha2[2,t] + 2*AA*CC^2*alpha1[4,t]*alpha2[2,t] - 4*AA*CC^2*alpha1[2,t]*alpha2[4,t] + 
        2*AA*CC^2*alpha1[2,t]*alpha3[2,t] - 4*AA*CC^2*alpha1[4,t]*alpha3[2,t] - 4*AA*CC^2*alpha1[4,t]*alpha3[4,t] + 6*AA*CC^2*alpha2[2,t]*alpha3[2,t] + 8*AA*CC^2*alpha2[2,t]*alpha3[4,t] + 
        8*AA*CC^2*alpha2[4,t]*alpha3[2,t] + 8*AA*CC^2*alpha2[4,t]*alpha3[4,t] - 16*BB*CC^2*alpha1[1,t]*alpha2[2,t] - 16*BB*CC^2*alpha1[3,t]*alpha2[2,t] + 8*BB*CC^2*alpha2[1,t]*alpha1[4,t] - 
        16*BB*CC^2*alpha1[1,t]*alpha2[4,t] - 16*BB*CC^2*alpha1[3,t]*alpha2[4,t] + 8*BB*CC^2*alpha2[3,t]*alpha1[4,t] - 16*BB*CC^2*alpha1[1,t]*alpha3[2,t] + 8*BB*CC^2*alpha3[1,t]*alpha1[2,t] - 
        16*BB*CC^2*alpha1[1,t]*alpha3[4,t] - 16*BB*CC^2*alpha1[3,t]*alpha3[2,t] + 8*BB*CC^2*alpha3[3,t]*alpha1[2,t] - 16*BB*CC^2*alpha1[3,t]*alpha3[4,t] + 3*CC^3*alpha1[1,t]*alpha2[2,t]*alpha3[2,t] + 
        CC^3*alpha2[1,t]*alpha1[2,t]*alpha3[2,t] - 2*CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[2,t] + 4*CC^3*alpha1[1,t]*alpha2[2,t]*alpha3[4,t] + 3*CC^3*alpha1[3,t]*alpha2[2,t]*alpha3[2,t] - 
        2*CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[2,t] + CC^3*alpha3[1,t]*alpha1[4,t]*alpha2[2,t] - 2*CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[2,t] + 4*CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[2,t] + 
        4*CC^3*alpha1[3,t]*alpha2[2,t]*alpha3[4,t] - 2*CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[4,t] + CC^3*alpha2[3,t]*alpha1[2,t]*alpha3[2,t] - 2*CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[4,t] + 
        CC^3*alpha3[3,t]*alpha1[4,t]*alpha2[2,t] + 4*CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[4,t] + 4*CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[2,t] - 2*CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[2,t] - 
        2*CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[4,t] + 4*CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[4,t] - 2*CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[4,t] + 16*AA*BB*CC*alpha1[2,t] + 16*AA*BB*CC*alpha1[4,t] - 
        32*AA*BB*CC*alpha2[2,t] - 32*AA*BB*CC*alpha2[4,t] - 32*AA*BB*CC*alpha3[2,t] - 32*AA*BB*CC*alpha3[4,t])) / (128*BB^2*CC*alpha1[2,t] - 512*BB^3 + 128*BB^2*CC*alpha1[4,t] + 128*BB^2*CC*alpha2[2,t] + 
        128*BB^2*CC*alpha2[4,t] + 128*BB^2*CC*alpha3[2,t] + 128*BB^2*CC*alpha3[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[4,t] - 
        24*BB*CC^2*alpha1[4,t]*alpha2[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha3[2,t] - 24*BB*CC^2*alpha1[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[4,t] - 
        24*BB*CC^2*alpha2[2,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[2,t] + 
        6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[2,t] + 9*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[2,t] + 9*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[4,t] + 
        6*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[2,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[4,t])
     
    
        br[2,t] = -(2*(128*AA*BB^2 + 64*BB^2*CC*alpha2[1,t] + 64*BB^2*CC*alpha2[3,t] - 4*AA*CC^2*alpha1[2,t]*alpha2[2,t] - 4*AA*CC^2*alpha1[4,t]*alpha2[2,t] + 2*AA*CC^2*alpha1[2,t]*alpha2[4,t] + 
        8*AA*CC^2*alpha1[2,t]*alpha3[2,t] + 6*AA*CC^2*alpha1[2,t]*alpha3[4,t] + 8*AA*CC^2*alpha1[4,t]*alpha3[2,t] + 8*AA*CC^2*alpha1[4,t]*alpha3[4,t] + 2*AA*CC^2*alpha2[2,t]*alpha3[4,t] - 
        4*AA*CC^2*alpha2[4,t]*alpha3[2,t] - 4*AA*CC^2*alpha2[4,t]*alpha3[4,t] - 16*BB*CC^2*alpha2[1,t]*alpha1[2,t] - 16*BB*CC^2*alpha2[1,t]*alpha1[4,t] + 8*BB*CC^2*alpha1[1,t]*alpha2[4,t] - 
        16*BB*CC^2*alpha2[3,t]*alpha1[2,t] + 8*BB*CC^2*alpha1[3,t]*alpha2[4,t] - 16*BB*CC^2*alpha2[3,t]*alpha1[4,t] - 16*BB*CC^2*alpha2[1,t]*alpha3[2,t] + 8*BB*CC^2*alpha3[1,t]*alpha2[2,t] - 
        16*BB*CC^2*alpha2[1,t]*alpha3[4,t] + 8*BB*CC^2*alpha3[3,t]*alpha2[2,t] - 16*BB*CC^2*alpha2[3,t]*alpha3[2,t] - 16*BB*CC^2*alpha2[3,t]*alpha3[4,t] + 4*CC^3*alpha2[1,t]*alpha1[2,t]*alpha3[2,t] - 
        2*CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[2,t] + CC^3*alpha1[1,t]*alpha2[2,t]*alpha3[4,t] + 3*CC^3*alpha2[1,t]*alpha1[2,t]*alpha3[4,t] + 4*CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[2,t] - 
        2*CC^3*alpha3[1,t]*alpha1[4,t]*alpha2[2,t] - 2*CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[2,t] - 2*CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[2,t] + CC^3*alpha1[3,t]*alpha2[2,t]*alpha3[4,t] + 
        4*CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[4,t] + 4*CC^3*alpha2[3,t]*alpha1[2,t]*alpha3[2,t] + CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[4,t] - 2*CC^3*alpha3[3,t]*alpha1[4,t]*alpha2[2,t] - 
        2*CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[4,t] - 2*CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[2,t] + 3*CC^3*alpha2[3,t]*alpha1[2,t]*alpha3[4,t] + 4*CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[2,t] + 
        CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[4,t] - 2*CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[4,t] + 4*CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[4,t] - 32*AA*BB*CC*alpha1[2,t] - 32*AA*BB*CC*alpha1[4,t] + 
        16*AA*BB*CC*alpha2[2,t] + 16*AA*BB*CC*alpha2[4,t] - 32*AA*BB*CC*alpha3[2,t] - 32*AA*BB*CC*alpha3[4,t])) / (128*BB^2*CC*alpha1[2,t] - 512*BB^3 + 128*BB^2*CC*alpha1[4,t] + 128*BB^2*CC*alpha2[2,t] + 
        128*BB^2*CC*alpha2[4,t] + 128*BB^2*CC*alpha3[2,t] + 128*BB^2*CC*alpha3[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[4,t] - 
        24*BB*CC^2*alpha1[4,t]*alpha2[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha3[2,t] - 24*BB*CC^2*alpha1[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[4,t] - 
        24*BB*CC^2*alpha2[2,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[2,t] + 
        6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[2,t] + 9*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[2,t] + 9*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[4,t] + 
        6*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[2,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[4,t])
     
     
        br[3,t] = -(2*(128*AA*BB^2 + 64*BB^2*CC*alpha3[1,t] + 64*BB^2*CC*alpha3[3,t] + 8*AA*CC^2*alpha1[2,t]*alpha2[2,t] + 8*AA*CC^2*alpha1[4,t]*alpha2[2,t] + 8*AA*CC^2*alpha1[2,t]*alpha2[4,t] + 
        6*AA*CC^2*alpha1[4,t]*alpha2[4,t] - 4*AA*CC^2*alpha1[2,t]*alpha3[2,t] - 4*AA*CC^2*alpha1[4,t]*alpha3[2,t] + 2*AA*CC^2*alpha1[4,t]*alpha3[4,t] - 4*AA*CC^2*alpha2[2,t]*alpha3[4,t] + 
        2*AA*CC^2*alpha2[4,t]*alpha3[2,t] - 4*AA*CC^2*alpha2[4,t]*alpha3[4,t] - 16*BB*CC^2*alpha3[1,t]*alpha1[2,t] + 8*BB*CC^2*alpha1[1,t]*alpha3[4,t] - 16*BB*CC^2*alpha3[1,t]*alpha1[4,t] - 
        16*BB*CC^2*alpha3[3,t]*alpha1[2,t] + 8*BB*CC^2*alpha1[3,t]*alpha3[4,t] - 16*BB*CC^2*alpha3[3,t]*alpha1[4,t] + 8*BB*CC^2*alpha2[1,t]*alpha3[2,t] - 16*BB*CC^2*alpha3[1,t]*alpha2[2,t] - 
        16*BB*CC^2*alpha3[3,t]*alpha2[2,t] + 8*BB*CC^2*alpha2[3,t]*alpha3[2,t] - 16*BB*CC^2*alpha3[1,t]*alpha2[4,t] - 16*BB*CC^2*alpha3[3,t]*alpha2[4,t] - 2*CC^3*alpha2[1,t]*alpha1[2,t]*alpha3[2,t] + 
        4*CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[2,t] - 2*CC^3*alpha1[1,t]*alpha2[2,t]*alpha3[4,t] - 2*CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[2,t] + 4*CC^3*alpha3[1,t]*alpha1[4,t]*alpha2[2,t] + 
        4*CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[2,t] + CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[2,t] - 2*CC^3*alpha1[3,t]*alpha2[2,t]*alpha3[4,t] + CC^3*alpha2[1,t]*alpha1[4,t]*alpha3[4,t] - 
        2*CC^3*alpha2[3,t]*alpha1[2,t]*alpha3[2,t] + 4*CC^3*alpha3[1,t]*alpha1[2,t]*alpha2[4,t] + 4*CC^3*alpha3[3,t]*alpha1[4,t]*alpha2[2,t] - 2*CC^3*alpha1[1,t]*alpha2[4,t]*alpha3[4,t] +
         CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[2,t] - 2*CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[2,t] + 3*CC^3*alpha3[1,t]*alpha1[4,t]*alpha2[4,t] + 4*CC^3*alpha3[3,t]*alpha1[2,t]*alpha2[4,t] - 
         2*CC^3*alpha1[3,t]*alpha2[4,t]*alpha3[4,t] + CC^3*alpha2[3,t]*alpha1[4,t]*alpha3[4,t] + 3*CC^3*alpha3[3,t]*alpha1[4,t]*alpha2[4,t] - 32*AA*BB*CC*alpha1[2,t] - 32*AA*BB*CC*alpha1[4,t] - 
         32*AA*BB*CC*alpha2[2,t] - 32*AA*BB*CC*alpha2[4,t] + 16*AA*BB*CC*alpha3[2,t] + 16*AA*BB*CC*alpha3[4,t])) / (128*BB^2*CC*alpha1[2,t] - 512*BB^3 + 128*BB^2*CC*alpha1[4,t] + 128*BB^2*CC*alpha2[2,t] + 
         128*BB^2*CC*alpha2[4,t] + 128*BB^2*CC*alpha3[2,t] + 128*BB^2*CC*alpha3[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha2[2,t] - 32*BB*CC^2*alpha1[2,t]*alpha2[4,t] - 
         24*BB*CC^2*alpha1[4,t]*alpha2[4,t] - 32*BB*CC^2*alpha1[2,t]*alpha3[2,t] - 24*BB*CC^2*alpha1[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha1[4,t]*alpha3[4,t] - 
         24*BB*CC^2*alpha2[2,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[2,t]*alpha3[4,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[2,t] - 32*BB*CC^2*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[2,t] + 
         6*CC^3*alpha1[2,t]*alpha2[2,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[2,t] + 9*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[2,t] + 9*CC^3*alpha1[4,t]*alpha2[2,t]*alpha3[4,t] + 
         6*CC^3*alpha1[2,t]*alpha2[4,t]*alpha3[4,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[2,t] + 6*CC^3*alpha1[4,t]*alpha2[4,t]*alpha3[4,t])
        
        p[1,t] = br[1,t] + wn[1,t]
        p[2,t] = br[2,t] + wn[2,t]
        p[3,t] = br[3,t] + wn[3,t]

        pavg[:,t+1]=pavg[:,t] +lambda*( p[:,t] -pavg[:,t]);
        
        p12_avg[t] = (p[1,t] + p[2,t]) / 2
        p13_avg[t] = (p[1,t] + p[3,t]) / 2
        p23_avg[t] = (p[2,t] + p[3,t]) / 2
            
        alpha1[1:2,t+1]=alpha1[1:2,t]+lambda*inv(R12)*[1; p13_avg[1,t]]*(p[2,t]-alpha1[1,t]-alpha1[2,t]*p13_avg[1,t]);
        R12 = R12 +lambda*( [1  p13_avg[1,t]; p13_avg[1,t] p13_avg[1,t]*p13_avg[1,t]] -R12);

        upper= 1; 

        if alpha1[2,t+1] > upper
            alpha1[2,t+1] = upper;
            alpha1[1,t+1] = pavg[2,t+1]-((pavg[1,t+1]+pavg[3,t+1])/2);

            elseif alpha1[2,t+1] <0
                alpha1[2,t+1]=0;
                alpha1[1,t+1]=pavg[2,t+1];
        end
    
        if ~all(isfinite.(R12))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R12) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
    
        alpha1[3:4,t+1]=alpha1[3:4,t]+lambda*inv(R13)*[1; p12_avg[1,t]]*(p[3,t]-alpha1[3,t]-alpha1[4,t]*p12_avg[1,t]);
        R13 = R13 +lambda*( [1  p12_avg[1,t]; p12_avg[1,t] p12_avg[1,t]*p12_avg[1,t]] -R13);
        
        if alpha1[4,t+1] > upper
            alpha1[4,t+1] = upper;
            alpha1[3,t+1]=pavg[3,t+1]-((pavg[1,t+1]+pavg[2,t+1])/2);

            elseif alpha1[4,t+1] <0
                alpha1[4,t+1]=0;
                alpha1[3,t+1]=pavg[3,t+1];
        end
    
        if ~all(isfinite.(R13))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R13) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
        
        alpha2[1:2,t+1]=alpha2[1:2,t]+lambda*inv(R21)*[1; p23_avg[1,t]]*(p[1,t]-alpha2[1,t]-alpha2[2,t]*p23_avg[1,t]);
        R21 = R21 +lambda*( [1  p23_avg[1,t]; p23_avg[1,t] p23_avg[1,t]*p23_avg[1,t]] -R21);
        
        if alpha2[2,t+1] > upper
            alpha2[2,t+1] = upper;
            alpha2[1,t+1]=pavg[1,t+1]-((pavg[2,t+1]+pavg[3,t+1])/2);

            elseif alpha2[2,t+1] <0
                alpha2[2,t+1]=0;
                alpha2[1,t+1]=pavg[1,t+1];
        end

        if ~all(isfinite.(R21))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R21) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
    
        alpha2[3:4,t+1]=alpha2[3:4,t]+lambda*inv(R23)*[1; p12_avg[1,t]]*(p[3,t]-alpha2[3,t]-alpha2[4,t]*p12_avg[1,t]);
        R23 = R23 +lambda*( [1  p12_avg[1,t]; p12_avg[1,t] p12_avg[1,t]*p12_avg[1,t]] -R23);
        
        if alpha2[4,t+1] > upper
            alpha2[4,t+1] = upper;
            alpha2[3,t+1]=pavg[3,t+1]-((pavg[1,t+1]+pavg[2,t+1])/2);

            elseif alpha2[4,t+1] <0
                alpha2[4,t+1]=0;
                alpha2[3,t+1]=pavg[3,t+1];
        end

        if ~all(isfinite.(R23))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R23) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
    
        alpha3[1:2,t+1]=alpha3[1:2,t]+lambda*inv(R31)*[1; p23_avg[1,t]]*(p[1,t]-alpha3[1,t]-alpha3[2,t]*p23_avg[1,t]);
        R31 = R31 +lambda*( [1  p23_avg[1,t]; p23_avg[1,t] p23_avg[1,t]*p23_avg[1,t]] -R31);

        if alpha3[2,t+1] > upper
            alpha3[2,t+1] = upper;
            alpha3[1,t+1]=pavg[1,t+1]-((pavg[2,t+1]+pavg[3,t+1])/2);

            elseif alpha3[2,t+1] <0
                alpha3[2,t+1]=0;
                alpha3[1,t+1]=pavg[1,t+1];
        end
    
        if ~all(isfinite.(R31))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R31) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
    
        alpha3[3:4,t+1]=alpha3[3:4,t]+lambda*inv(R32)*[1; p13_avg[1,t]]*(p[2,t]-alpha3[3,t]-alpha3[4,t]*p13_avg[1,t]);
        R32 = R32 +lambda*( [1  p13_avg[1,t]; p13_avg[1,t] p13_avg[1,t]*p13_avg[1,t]] -R32);
    
        if alpha3[4,t+1] > upper
            alpha3[4,t+1] = upper;
            alpha3[3,t+1]=pavg[2,t+1]-((pavg[1,t+1]+pavg[3,t+1])/2);

            elseif alpha3[4,t+1] <0
                alpha3[4,t+1]=0;
                alpha3[3,t+1]=pavg[2,t+1];
        end

        if ~all(isfinite.(R32))
            #warning = ('Matrix is Singular');
            break
        else
            if rank(R32) < 2
                #warning = ('Matrix is Singular');
                break
            end
        end
        
    end
    return vcat(p,alpha1,alpha2,alpha3);
end