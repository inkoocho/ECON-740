clear;
close all;
rand('seed',sum(100*clock));
Toler=1e-6;
N0=0;
Nsim=20;
Tsim=3000;
Titer=Tsim + N0;
alpha=0.99;
rho=.95;
if (alpha*rho) > .933
  absc = 1 - 16*alpha*rho*(1-alpha*rho);
  pibar = (4*alpha*rho - 1 + sqrt(absc))/(4*alpha*rho);
  pibar
end
sigmap=0.2; sigmab=0.0001;
%sigmaf(1) = .000005; sigmaf(2) = .0001; sigmaf(3) = .001; sigmaf(4)=0.01;
%Numf = length(sigmaf);

sigmafmin=0.0000025;   sigmafmax=0.0005; sigmafinc=0.1*(sigmafmax-sigmafmin);
sigmaf=sigmafmin: sigmafinc: sigmafmax ;
Numf=length(sigmaf);

pihatgain=0.00000001;
delta = (1-rho*alpha);
gam = (1-alpha);
sce = delta/(1-alpha*rho);
result =zeros(Numf,Nsim);

for k=1:Numf
Vbar=sqrt(sigmap*sigmab*(1-rho^2))/sqrt(sigmaf(k));
gain=sqrt(sigmab*sigmaf(k)/(sigmap*(1-rho^2)));
V0(1) = Vbar;
V1(1) = Vbar;
convg = 0;
for j=1:Nsim
  price(1)=gam/(1-alpha);
  beta0(1) = gam/(1-alpha) + .01*randn*sqrt(V0(1));
  beta1(1) = price(1) + .01*randn*sqrt(V1(1));
  pi(1) = 0.5; r = 0.0; pihat(1) = 0.5;
  %pi(1) = 0.0; r = 0.0; pihat(1) = 0.0; 

  vp=sqrt(sigmap)*randn(Titer);
  vf=sqrt(sigmaf(k))*randn(Titer);
  g=sqrt(sigmab*sigmaf(k))/sqrt(sigmap*(1-rho^2));
  f(1) = vf(1);
  A1(1)=1; A0(1)=1;
  %phi(:,1) = [1; f(1)];
  beta0v(:,1) = [gam/(1-alpha); beta0(1)];
  sumvarp = 0.0;

  for i=2:Titer
      
      f(i) = rho*f(i-1) + vf(i);

      res1 = price(i-1) - beta1(i-1);   
  %   res0 = price(i-1) - (pihat(i-1)*beta1(i-1)+(1-pi(i-1))*beta0(i-1)) - sce*f(i-1); %(7.22)
  %   res0 = price(i-1) - (pihat(i-1)*beta1(i-1)+(1-pi(i-1))*(beta0(i-1) + sce*f(i-1))); 
      res0 = price(i-1) - (pihat(i-1)*beta1(i-1)+(1-pihat(i-1))*(beta0(i-1) + sce*f(i-1)));

      beta1(i) = beta1(i-1) + (V1(i-1)/(sigmap + V1(i-1)))*res1;

      V1(i) = V1(i-1) - ((V1(i-1))^2)/(sigmap + V1(i-1)) + sigmab;
      V0(i) = V0(i-1) - ((V0(i-1))^2)/(sigmap + V0(i-1));

      A1(i) = (exp(-.5*(res1^2)/(sigmap + V1(i-1))))/sqrt(2*3.14159*(sigmap + V1(i-1)));
      
      M0gaincp = V0(i-1)/(sigmap + V0(i-1));
   
      beta0(i) = beta0(i-1) + M0gaincp*res0;
      
      A0(i) = (exp(-.5*(res0^2)/(sigmap + V0(i-1))))/sqrt(2*3.14159*(sigmap + V0(i-1)));
      
      r(i) = r(i-1) + log(A1(i)/A0(i));
      pi(i) = 1/(1 + exp(-r(i)));
      pihat(i) = pihat(i-1) + (1/(50+i))*(pi(i) - pihat(i-1));
      %pihat(i) = 0.0;
      
      price(i) = gam + delta*f(i) + alpha*((pi(i) + (1-pi(i))*pihat(i))*beta1(i) + ...
          (1-pi(i))*(1-pihat(i))*beta0(i) + (1-pi(i))*(1-pihat(i))*rho*sce*f(i)) + vp(i);
  end
  if 1-pi(Titer) < pi(Titer)
      convg = convg+1;
  else
      convg = convg;
  end
  result(k,j)=pi(Titer);
  varp1 = var(price);
  sumvarp = sumvarp + varp1;
end
prop(k) = convg/Nsim;
varp(k) = sumvarp/Nsim;
stdevp(k) = sqrt(varp(k));
%  result(k,Nsim)=pi(Titer);
end

%%  End simulations

%figure;
%plot(sigmaf,prop); title('TVP Convergence Probability');

for ii=1:Numf
    av(ii)=mean(result(ii,:));
end

figure;
plot(sigmaf,av); title('TVP Convergence Probability');


xax=1:Tsim;
%beta0v1 = beta0v(2,1:Tsim);
figure;
subplot(2,2,1); plot(xax,price); title('Price');
subplot(2,2,2); plot(xax,pi); title('Pi'); 
subplot(2,2,3); plot(xax,beta0); title('Beta0');
subplot(2,2,4); plot(xax,beta1); title('Beta1');
figure;
plot(xax,pihat); title('Pihat');
stdevp

print -depsc fig1.ps





