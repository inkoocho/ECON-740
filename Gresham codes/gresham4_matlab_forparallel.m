function [pi,beta0,beta1,price]=gresham4_matlab_forparallel()
%clear;
%close all;
%rand('seed',sum(100*clock));
%global s
%rng(s)
%randn('seed',12345678);
Toler=1e-6;
N0=0;
Nsim=4000;
Titer=Nsim + N0;
alpha=0.96;
rho=.95;
sigmap=1.0; sigmaf=1.0; sigmab=0.0005;
delta = (1-rho*alpha);
sce = delta/(1-alpha*rho);

pi=ones(Titer,1);
beta0=ones(Titer,1);
beta1=ones(Titer,1);
price=ones(Titer,1);

V0(1) = 0.2;
V1(1) = 0.2;

beta0(1) = sce + randn*sqrt(V0(1));
beta1(1) = sce + randn*sqrt(V1(1));
pi(1) = 0.1; r = 0.0;


vp=sqrt(sigmap)*randn(Titer,1);
vf=sqrt(sigmaf)*randn(Titer,1);
g=sqrt(sigmab*sigmaf)/sqrt(sigmap*(1-rho^9));
varbeta = (g*sigmap)/(sigmaf*(1-alpha)*(2 - g*(1-alpha)));
varf = sigmaf/(1-rho^2);
varratio = 100*(varbeta*varf)/(varf*sce^2 + sigmap);
%delta
%sce
%g
%varbeta
%varratio
price(1)=0;
f(1) = 0.0;
A1(1)=1; A0(1)=1;
for i=2:Titer

f(i) = rho*f(i-1) + vf(i);

res1 = price(i-1) - beta1(i-1)*f(i-1);
res0 = price(i-1) - beta0(i-1)*f(i-1);
MSE1(i-1) = sigmap + V1(i-1)*f(i-1)^2 + f(i-1)*f(i-1)*V1(i-1)*(1-alpha*rho*pi(i-1))^2;

MSE0(i-1) = sigmap + ((alpha*rho*pi(i-1))^2)*V1(i-1)*f(i-1)^2 + (f(i-1)^2)*V0(i-1)*(1-alpha*rho*(1-pi(i-1)))^2+V0(i-1)*(f(i-1)^2);%%% added the last term

beta1(i) = beta1(i-1) + (V1(i-1)/MSE1(i-1))*f(i-1)*res1;
beta0(i) = beta0(i-1) + (V0(i-1)/MSE0(i-1))*f(i-1)*res0;
V1(i) = V1(i-1) - ((f(i-1)*V1(i-1))^2)/MSE1(i-1) + sigmab; %%new Nov.7
V0(i) = V0(i-1) - ((f(i-1)*V0(i-1))^2)/MSE0(i-1) + 0; %%new Nov.7
%V1(i) = V1(i-1) - ((f(i-1)*V1(i-1))^2)/(sigmap + V1(i-1)*f(i-1)^2) + sigmab;
%V0(i) = V0(i-1) - ((f(i-1)*V0(i-1))^2)/(sigmap + V0(i-1)*f(i-1)^2) + 0;
A1(i) = exp(-.5*(res1^2)/MSE1(i-1))/sqrt(2*3.14159*MSE1(i-1));
A0(i) = exp(-.5*(res0^2)/MSE0(i-1))/sqrt(2*3.14159*MSE0(i-1));
r(i) = r(i-1) + log(A1(i)/A0(i));
pi(i) = 1/(1 + exp(-r(i)));
price(i) = (delta + alpha*rho*(pi(i)*beta1(i) + (1-pi(i))*beta0(i)))*f(i) + vp(i);

end

gresham_path=[pi,beta0,beta1,price];
end






