%%%% 01.09.2014 V2.1 Beta for: SVHJ; Alpha for: SVJ 
%% Intro
clc;clear;

%% Simulation parameters
  muj0=-0.05;
  mujq0=-0.14;
  sigmaj0=0.06;
  eta0=2.4;
  k0=4.8;
  vbar0=0.01;
  sigma0=0.22;
  rho0=-0.6;
  kl0=15; %Starting value is the difference from delta
  lambda0=5;
  delta0=5;
  theta = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0];%theta SVHj
% theta = [eta0,k0,vbar0,sigma0,rho0]; %theta SV
% theta_SVJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,lambda0]; %theta SV

%% Simulate data-set of states & derivatives
  seed=1+floor(1e6*rand(1,1));
  N=250;
  states=sim_states([log(100),0.2^2,5],N,1/250,'SVHJ',theta,{[100],[seed]});
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   tic
%   [deriv_data, true_states] = sim_deriv_data([log(100),0.2^2,5],N,5/250,'SV',theta,[100],[seed]);
%   toc
%   deriv_data(:,7)=deriv_data(:,5);
%   fprintf('\n')
  %Determine the structure of the data file
%         M=max(deriv_data(:,1));
%         index=1;
%         i=1;    
%         in_day=0;
%         while i<=M
%           TT_marks(i,1)=i;
%           TT_marks(i,2)=0;
%           %Pick the options in each day
%           in_day=1;
%           in_TT=3;
%           TT_now=deriv_data(index,2);
%           while index<=length(deriv_data) && deriv_data(index,1)<=i
%           if(TT_now<deriv_data(index,2))
%               TT_marks(i,in_TT)=in_day-1;
%               in_TT=in_TT+1;
%               TT_now=deriv_data(index,2);
%           end
%           index=index+1;
%           in_day=in_day+1;
%           end  
%           TT_marks(i,in_TT)=in_day-1;
%           i=i+1;
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end_OLD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate parameters based on the simulated sample
%     lb=[0, -0.40, 0.02, 0.5,   2, 0.005, 0.1, -0.8, 0.5, 0.05, 3]; %bound for theta(1) and theta(9) are here bounds for the difference
%     ub=[0.1, -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3, 30, 1.50, 30];
%     isnotokay=1;
%     fprintf('Picking theta_start\n');
%     while isnotokay
%     thetas=theta+((-1).^round(rand(1,11))).*0.2.*(1-rand(1,11)).*theta;
%     if sum(imag(boundtofree('SVHJ',thetas,lb,ub)))==0
%         isnotokay=0;
%     end
%     end
%     fprintf('Theta was picked\n');
%     theta=thetas;
%     sim_result=estimate2( deriv_data, TT_marks, 5/250,'SVHJ', theta,lb,ub,theta);
%     name_for_save=sprintf('deriv_data_%3.0f_%5.0f.mat',N,seed);
%     save(sprintf('Data Sim/Deriv Samples/%s',name_for_save));
    
    
    
    lb=[0, -0.40, 0.02, 0.5,   2, 0.01, 0.1, -0.8,0.1]; %lb SV
    ub=[0.07, -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3,5 ]; %ub SV
% theta = lb+rand([1,11]).*(ub-lb);
% theta(9) = theta(11) + theta(9); %theta(9) stores the difference, theta(11) the level
% theta(1) = theta(2) + theta(1); %theta(1) stores the difference, theta(2) the level
test(:,1) = 0.8:0.1:2;
bad_start=theta;%+0.2.*(1-rand(1,5)).*theta ;   %%%+[0.6.*(1-rand(1,3)).*theta(1:3),zeros([1,8])];
for q=1:length(test(:,1))
mod=bad_start;
mod(4)=theta(4)*test(q,1);
test(q,2)=quickcrit_true_states(states,5/250,'SV', mod,lb,ub,bad_start);
end
plot(test(:,1),test(:,2));

% % % deriv_data(:,4)=deriv_data(:,4)+normrnd(0,0.001,size(deriv_data(:,4)));

% % 
% theta1step=theta;
% X=impl_state_func_fut(deriv_data,TT_marks,'SVHJ',theta);
% plot([X(:,3),true_states(:,3)])
%   factor=1e3;
%   alpha=0.01;
% H0=funcH0(theta1step, X, 5/250, 'SVHJ', factor);
% fun = @(theta) estimate_hess(deriv_data, TT_marks, 5/250,'SVHJ',theta',lb,ub,theta1step,X,H0);
% out=hessian_2sided(fun,theta');