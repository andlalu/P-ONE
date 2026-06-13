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
delta0 = 10;

theta_SV = [eta0,k0,vbar0,sigma0,rho0]; %theta SV
theta_SVJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,lambda0]; %theta SVJ
theta_SVHJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ

%% Simulate data-set of states & derivatives
%  Set-up the simulated sample size  
seed = 1+floor(1e6*rand(1,1));
N = 501; %Sample size simulated length
dt=5/250;
type='SVHJ';
%Simulate states and derivatives prices:
[deriv_data, states, TT_marks, coeff_mat ] = sim_deriv_data([log(100),0.2^2,5],N,dt,type,theta_SVHJ,1111);%640449

%% Set-up routine characteristics
c=100;
theta_SVHJ2=theta_SVHJ;
theta_SVHJ2(10)=theta_SVHJ(10)./c;
theta_SVHJ2(11)=theta_SVHJ(11)./c;  
theta_SVHJ3=theta_SVHJ2;
% % % % % % % % % % % % % % theta_SVHJ3(9)=rr23;%!!!!!!!!!!!!!!%
lb=[-.4, -0.40, 0.02, 0.5,   2, 0.005, 0.1, -0.8, 0.5, 0.05./c, 3./c]; %bound for theta(1) and theta(9) are here bounds for the difference
ub=[0,  -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3, 30, 1.50./c, 30./c];
% states_imp_aux=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ3.*[ones(1,9) c c]).*repmat([1 1 1./c],[N-1,1]);
aux=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ3.*[ones(1,9) c c]).*repmat([1 1 1./c],[N-1,1]);
[OptMatrix, H0]=funcH0(theta_SVHJ3,aux, dt, type, 1e-18);
tic
%% Parameter-by-parameter plot 

out1=[]
out2=[]
ratios=0.6:.2:1.4;
states_imp=zeros([N-1,3,length(ratios)]);
vartest=7;
for rate=ratios
    rate
theta_SVHJ3(vartest)=theta_SVHJ2(vartest).*rate;
states_imp(:,:,ratios==rate)=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ3.*[ones(1,9) c c]).*repmat([1 1 1./c],[N-1,1]);
% % out1=[out1; critstep11(theta_SVHJ3, states_imp(:,:,ratios==rate), dt, type)];
out2=[out2; critstep2(theta_SVHJ3,  states_imp(:,:,ratios==rate), dt, type, H0, OptMatrix)];
end
aux=1.*(ratios==1);aux(aux==0)=NaN;
% subplot(2,2,1)
% plot(ratios.*theta_SVHJ(vartest),out1,ratios.*theta_SVHJ(vartest),aux.*min(out1),'*k')
% legend('1st step')
subplot(2,2,2)
plot(ratios.*theta_SVHJ(vartest),out2,'r',ratios.*theta_SVHJ(vartest),aux.*min(out2),'*k')
legend('2nd step')
subplot(2,2,3)
plot(squeeze(states_imp(:,2,:))-repmat(states(2:end,2),[1,length(ratios)]))
legend('Implied vols')
subplot(2,2,4)
plot(squeeze(states_imp(:,3,:))-repmat(states(2:end,3)./c,[1,length(ratios)]))
legend('Implied intens')


%% Minimization routine
% [states_imp,x]=impl_state_func_fut(deriv_data, TT_marks, type,[theta_SVHJ2(1:9) theta_SVHJ2(10)*c theta_SVHJ2(11)*c]);
% subplot(1,2,1)
% plot([states(2:end,2)-states_imp(:,2),states(2:end,2)]);
% subplot(1,2,2)
% plot([states(2:end,3)-states_imp(:,3),states(2:end,3)]);
% states_estim = [diff(states(:,1)),states_imp2(2:end,2),states_imp2(2:end,3)./c]; %create stationary return series for estimation
% [OptMatrix, H0]=funcH0(theta_SVHJ2, states_estim, dt, type, alpha);
% func2 = @(theta) critstep1(freetobound(type,theta,lb,ub),states_estim,dt,type);%,H0,OptMatrix);
% % % func3 = @(theta_free) func2(freetobound(type,theta_free,lb,ub));
% func2(theta_SVHJ)
% options=optimset('Display','iter','TolFun',1e-2,'MaxIter',150);
% aux=fminsearch(func3,boundtofree(type,theta_SVHJ2,lb,ub),options)
% out2(i,:)=freetobound(type,aux,lb,ub)
% aux=fminsearch(func2,boundtofree(type,theta_SVHJ2,lb,ub),options)
% out2(i,:)=freetobound(type,aux,lb,ub)
% norm2=max(abs((theta_SVHJ2-out2(i,:))))
% theta_SVHJ2=out2(i,:);
% i=i+1
toc
%% Make some noise at end of routine
load('chirp');
% load('handel');
sound(y);
disp('Done!');  