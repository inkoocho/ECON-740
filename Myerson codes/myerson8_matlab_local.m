
% produce an estimation of the probability of deviating from optimal price
% by epsilon, using an adaptive algm. 
% Originally written by Cheng Ding in 2021.
%
% given parameters, simulate

clear
clc
close all
% rand('seed',sum(100*clock));
tic;

mu = 10; % mean of the trucated normal dist of consumer's value.  Fixed.
% sigma = 16; % std of the above dist 9 is minimum.   
sigmamin =11; % if less than 9, everyone buys.
sigmaincre =0.001;  % increment of sigma
jsigma=5000;   % number of demand curves.   Generate different demands by changing sigma.

Nbuyer = 100; % # of buyers in each period 
a = 0.0001; % learning rate
epsilon = 0.75; % size of perturbation to pricing
T= 300000; % # of experiments for each demand curve

N=T;   % tentative    Need to be cleaned up.

CRsample=5;

% Tlist = 500:500:50000; 

tol = 4e-1; % accuracy bound for the fzero to calculate the true optimal price

% store the value

truemu =zeros(1,jsigma);
truesigma = zeros(1,jsigma);
truepstar = zeros(1,jsigma);
truepistar = zeros(1,jsigma);
estimated =zeros(2,T+1,jsigma);
forecastp =zeros(T,jsigma);
realizedp =zeros(T,jsigma);
realizedq =zeros(T,jsigma);

CLferror=zeros(1,jsigma);
CRferror=zeros(CRsample,jsigma);
estimatedcr=zeros(CRsample,jsigma);

% number of demand curves which generates a good prediction after T rounds.
count = 0;

parfor j = 1:jsigma     % Generate different demand curve by changing sigma

    sigma =sigmamin + sigmaincre*(j-1);

    hazardrate = @(x) hazard(x, mu, sigma, mu) - x;
    [pstar fval] = fzero(hazardrate, [mu, 3*mu]);
    if abs(fval) > 1e-3
        warning('optimal pricing is not found given current parametrization')
    end
    pistar = pstar*(1 - truncated_normal_a_cdf(pstar, mu, sigma, mu));

    truemu(j)=mu;
    truesigma(j)=sigma;
    truepstar(j)=pstar;
    truepistar(j)=pistar;

%pmat = zeros(N,length(Tlist));
%Prob = zeros(1,length(Tlist));

%f = waitbar(0, 'Starting');
%parfor i = 1:length(Tlist)
%  i=1;
% T = Tlist(i);
%    display(i)
    
            beta = zeros(2, T);
            beta(1,1) = 15; 
            beta(2,1)= -0.5;   % initial estimated price is beta0
                        
            R = zeros(2, 2, T);
            R(:,:,1) = eye(2);

            rand('seed',sum(100*clock));   % initialize the random number generator for each demand

            e = 2*epsilon*(binornd(1, 0.5, [1,T]) - 0.5);
            eprice=zeros(1,T);
            p = zeros(1, T);
            Fp = zeros(1, T);
            Dagg = zeros(1, T);
            q = zeros(1, T);
            pi = zeros(1, T);
            
                      
                      
            beta(1,1) = 10; 
            beta(2,1)= -5;   %
            
            for iter = 1:T
                % pricing decision based on belief from last period
                eprice(iter)=-beta(1,iter)/(2*beta(2,iter));
                p(iter) =eprice(iter) + e(iter);
                
                % demand given pricing decision
                Fp(iter) = truncated_normal_a_cdf(p(iter), mu, sigma, mu);
                Dagg(iter) = binornd(Nbuyer, 1 - Fp(iter));
                q(iter) = Dagg(iter)/Nbuyer;
                
                %record the market outcome
                
                           
                pi(iter) = q(iter)*p(iter);
           
                if q(iter)==0
                    R(:,:,iter) = [1,eprice(iter);eprice(iter),eprice(iter)*eprice(iter)+epsilon*epsilon];
                    beta(:,iter+1) = beta(:,iter) + a*([-epsilon; 0] - beta(:,iter));
%                   R(:,:,iter+1) = R(:,:,iter) + a*([1,p(iter);p(iter),p(iter)^2] - R(:,:,iter));

                elseif q(iter)==1  
                    R(:,:,iter) = [1,eprice(iter);eprice(iter),eprice(iter)*eprice(iter)+epsilon*epsilon];                   
                    beta(:,iter+1) = beta(:,iter) + a*([epsilon; 0] - beta(:,iter));
%                   R(:,:,iter+1) = R(:,:,iter) + a*([1,p(iter);p(iter),p(iter)^2] - R(:,:,iter));
                else
                    R(:,:,iter) = [1,eprice(iter);eprice(iter),eprice(iter)*eprice(iter)+epsilon*epsilon];  
                    beta(:,iter+1) = beta(:,iter) + a*(R(:,:,iter)\([1;p(iter)]*(q(iter) - [1,p(iter)]*beta(:,iter))));
%                   R(:,:,iter+1) = R(:,:,iter) + a*([1,p(iter);p(iter),p(iter)^2] - R(:,:,iter));
                end
                
            end
 
            
 %       pmat(j,i) = p(T);   
%        if abs(p(T)-pstar) < tol
        % if abs( CLferror(j)) < tol
        %    count = count + 1;
        %    fprintf('good forecast error is %f \n',CLferror(j));
        %else
        %    fprintf('bad forecast error is %f \n',CLferror(j));
        %end            
            
 % record for the later use       
             
 CLferror(1,j)=eprice(T)-pstar;
 forecastp(:,j)=eprice;
 realizedp(:,j)=p;
 realizedq(:,j)=q;
                  
 estimated(:,:,j)=beta;
 %estimated(2,T,j)=beta(2,T);
 
 %%%Cole and Roughgarden  
    for isample=1:CRsample
        Obvalue=2*T*isample;
        
  value = truncated_normal_a_sample(mu,sigma,mu,Obvalue);
        q = (1:Obvalue)/Obvalue;
        pcr = sort(value, 'descend');
        % revenue given decreasing price sequence
        rev = q.*pcr;
        [M, I] = max(rev);
 CRferror(isample,j)= pstar -pcr(I);
 estimatedcr(isample,j)=pcr(I);
    end
        
end
    
toc   
 
%mean(CLferror)
%var(CLferror)
%hist(CLferror)

%mean(CRferror)
%var(CRferror)
%hist(CRferror)

%figure
%histogram(CLferror,-1.0:.01:1.0)
%hold on
%for isample=1:CRsample
%histogram(CRferror(isample,:),-1.0:.01:1.0)
%end
%box off
%axis tight
%legalpha('CL','CR','location','northwest')
%legend boxoff






%    end
% plot(p);
% plot(pi);

save("myerson8.mat", '-v7.3');
