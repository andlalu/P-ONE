%% Compare the Runge-Kutta solution of the PAN model with the solution in the paper for the CCF under the P measure:
clc; clear; 

%% Parameter set-up 
theta=[ -0.0387 -0.217  0.037    2.89   4.22    0.014    0.345   -0.45   28.12 ];

muj=theta(1);
mujq=theta(2);
sigmaj=theta(3);
eta=theta(4);
k=theta(5);
vbar=theta(6);
sigma=theta(7);
rho=theta(8);
lambda1=theta(9);
lambda0=0;
%% CCF value vector evaluation 
u_y=1;
u_v=1;
v=0.1^2;
Delta=5/250; %time period

%% First program the PAN paper solution: 
a = -u_y^2-2*u_y*(eta-0.5-lambda1*(exp(mujq+0.5*sigmaj^2)-1))-2*lambda1*(exp(u_y*muj+0.5*u_y^2*sigmaj^2)-1);

b = sigma*rho*u_y-k;

gamma = sqrt(b^2+a*sigma^2);

B = - (a*(1-exp(-gamma*Delta)) - u_v*(2*gamma-(gamma-b)*(1-exp(-gamma*Delta))))/(2*gamma-(gamma+b)*(1-exp(-gamma*Delta))-u_v*sigma^2*(1-exp(-gamma*Delta)))

A = -(k*vbar)/(sigma^2) *((gamma+b)*Delta + 2*log(1- (gamma+b+sigma^2*u_v)/(2*gamma)*(1-exp(-gamma*Delta))))+ (exp(u_y*muj+0.5*u_y^2*sigmaj^2)-1-u_y*(exp(mujq+0.5*sigmaj^2)-1))*lambda0*Delta

CCF_anlt=exp(A+B*v)

%% Compile the matrices in the DPS form for the Runge-Kutta evaluation part:
%Define CCF under P system matrices
K0 = [ -(exp(mujq+0.5*sigmaj^2)-1)*lambda0;      k*vbar;];
K1 = [ 0, eta-0.5-(exp(mujq+0.5*sigmaj^2)-1)*lambda1;
    0,  -k;];
H0 = [ 0,      0;
    0,      0;];
H1(:,:,1) = [ 0,0; 0,0;];
H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
L0=lambda0;
L1=[0;1;];
JT=@(beta) (exp(beta(1)*muj+0.5*beta(1)^2*sigmaj^2)-1);

coef_num = ccf(-1i.*[u_y;u_v],Delta,K0,K1,H0,H1,L0,L1,JT)

CCF_num= exp(sum([1,0,v].*coef_num))
