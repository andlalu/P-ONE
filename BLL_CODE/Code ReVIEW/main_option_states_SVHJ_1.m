%%%% 01.06.2015 V3.0 Beta for SVHJ
%% Intro
tic
clc;
clear;
%% Simulation parameters
muj0 = -0.05;
mujq0 = -0.14;
sigmaj0 = 0.06;
eta0 = 2.4;
k0 = 4.8;
vbar0 = 0.01;
sigma0 = 0.22;
rho0 = -0.6;
kl0 = 18; %Starting value is the difference from delta
lambda0 = 0.3;
delta0 = 10.5;

theta_SV = [eta0,k0,vbar0,sigma0,rho0]; %theta SV
theta_SVJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,lambda0]; %theta SVJ
theta_SVHJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ

%% Simulate data-set of states & derivatives
%  Set-up the simulated sample size  
seed = 1+floor(1e6*rand(1,1));
N = 501; %Sample size simulated length

% % %  SV Case - Simulate the state vector
dt=5/250;
type='SVHJ';

%%Parameter by parameter plots
%Simulate states and derivatives prices:
[deriv_data, states, TT_marks, coeff_mat ] = sim_deriv_data([log(100),0.2^2,5],N,dt,type,theta_SVHJ,848496);%640449
%Imply states back from prices assuming correct parameter vector:

%Temporary checking out criterion function as a function of one parameter

c=10;
alpha=0.02;

theta_SVHJ2=theta_SVHJ;
theta_SVHJ2(10)=theta_SVHJ(10)./c;
theta_SVHJ2(11)=theta_SVHJ(11)./c;
theta_aux=theta_SVHJ2;
theta_SVHJ2(9)=26;
theta_SVHJ2(11)=22./c;
lb=[-.4, -0.40, 0.02, 0.5,   2, 0.005, 0.1, -0.8, 0.5, 0.05./c, 3./c]; %bound for theta(1) and theta(9) are here bounds for the difference
ub=[0,  -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3, 30, 1.50./c, 30./c];

tic

func2 = @(theta) critstep1(theta,impl_state_func_fut(deriv_data, TT_marks, type,[theta(1:9) theta(10)*c theta(11)*c]).*repmat([1 1 1./c],[N-1,1]),dt,type);%,H0,OptMatrix);
func3 = @(theta_free) func2(freetobound(type,theta_free,lb,ub));
func2(theta_aux)
options=optimset('Display','iter','TolFun',1e-2,'MaxIter',400);
minimized=fminunc(func3,boundtofree(type,theta_SVHJ2,lb,ub),options)
freetobound(type,minimized,lb,ub)

toc
load('handel');
sound(y);
disp('Done!');  

%%%Phased OUT
% % [states_imp2,x]=impl_state_func_fut(deriv_data, TT_marks, type,[theta_SVHJ2(1:9) theta_SVHJ2(10)*c theta_SVHJ2(11)*c]);
% % states_estim = [diff(states(:,1)),states_imp2(2:end,2),states_imp2(2:end,3)./c]; %create stationary return series for estimation
% % [OptMatrix, H0]=funcH0(theta_SVHJ2, states_estim, dt, type, alpha);
% % out2(i,:)=freetobound(type,aux,lb,ub)
% % aux=fminsearch(func2,boundtofree(type,theta_SVHJ2,lb,ub),options)
% % out2(i,:)=freetobound(type,aux,lb,ub)
% % norm2=max(abs((theta_SVHJ2-out2(i,:))))
% % theta_SVHJ2=out2(i,:);
% % i=i+1