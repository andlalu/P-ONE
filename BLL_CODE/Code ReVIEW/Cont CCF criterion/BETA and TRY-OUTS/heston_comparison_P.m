%% Compare the Runge-Kutta solution of the SV model with the solution in the AFFINE BOOK MANUSCRIPT for the CCF under the P measure:
clc; clear; 

%% Parameter set-up 
  eta=2.4;
  k=4.8;
  vbar=0.01;
  sigma=0.22;
  rho=-0.6;
%% CCF value vector evaluation 
u_y=1;
u_v=1;
v=0.1^2;
Delta=5/250; %time period

%% First program the PAN paper solution: 
a = -u_y^2-2*u_y*(eta-0.5);

b = sigma*rho*u_y-k;

gamma = sqrt(b^2+a*sigma^2);

B = - (a*(1-exp(-gamma*Delta)) - u_v*(2*gamma-(gamma-b)*(1-exp(-gamma*Delta))))/(2*gamma-(gamma+b)*(1-exp(-gamma*Delta))-u_v*sigma^2*(1-exp(-gamma*Delta)))

A = -(k*vbar)/(sigma^2) *((gamma+b)*Delta + 2*log(1- (gamma+b+sigma^2*u_v)/(2*gamma)*(1-exp(-gamma*Delta))));

CCF_anlt=exp(A+B*v)


%% Compile the matrices in the DPS form for the Runge-Kutta evaluation part:
%Define CCF under P system matrices
K0 = [ 0;      k*vbar;];
K1 = [ 0, eta-0.5;
    0,  -k;];
H0 = [ 0,      0;
    0,      0;];
H1(:,:,1) = [ 0,0; 0,0;];
H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
L0=0;
L1=[0;0;];
JT=@(beta) (0);

coef_num = ccf(1i.*[u_y;u_v],Delta,K0,K1,H0,H1,L0,L1,JT)

CCF_num= exp(sum([1,0,v].*coef_num))
