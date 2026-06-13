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
delta0 = 16.5;

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
[deriv_data, states, TT_marks, coeff_mat ] = sim_deriv_data([log(100),0.2^2,5],N,dt,type,theta_SVHJ,1);%640449
%Imply states back from prices assuming correct parameter vector:

comp1=[]
comp2=[]
comp3=[]
% theta_SVHJ2=theta_SVHJ;%[3.14328830986134,6.28657661972268,0.0130970346244222,0.288134761737289,-0.785822077465335];%theta_SV*(0.5+rand(1,1));
theta_SVHJ2=[-0.1   -0.2    0.1    1    6    0.015    0.15   -0.5   25   0.3   10];
ratios=0.6:.1:1.4;
vartest=11
% states_imp=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ2);
% states_estim = [diff(states(:,1)),states_imp(2:end,2:end)]; %create stationary return series for estimation
alpha=1e-20;
states_imp=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ2);
states_estim = [diff(states(:,1)),states_imp(2:end,2:end)]; %create stationary return series for estimation
[OptMatrix1, H01]=funcH0(theta_SVHJ2,states_estim, dt, type, alpha);
i=1;
for rate=ratios
    theta_SVHJ2(vartest)=theta_SVHJ(vartest).*rate;
    rate
        states_imp=impl_state_func_fut(deriv_data, TT_marks, type, theta_SVHJ2);
        states_estim(:,:,i) = [diff(states(:,1)),states_imp(2:end,2:end)]; %create stationary return series for estimation
        [OptMatrix, H0]=funcH0(theta_SVHJ2,states_estim, dt, type, alpha);
        comp1 = [comp1, critstep1(theta_SVHJ2, states_estim, dt, type)];
        comp2 = [comp2, critstep2(theta_SVHJ2, states_estim, dt, type, H0, OptMatrix)];
        comp3 = [comp3, critstep2(theta_SVHJ2, states_estim, dt, type, H01, OptMatrix1)];
        i=i+1;
end
aux=1.*(ratios==1);aux(aux==0)=NaN;
subplot(1,3,1)
plot(ratios.*theta_SVHJ(vartest),comp1,ratios.*theta_SVHJ(vartest),aux.*min(comp1),'*k')
legend('1st step')
subplot(1,3,2)
plot(ratios.*theta_SVHJ(vartest),comp2,'r',ratios.*theta_SVHJ(vartest),aux.*min(comp2),'*k')
legend('2nd step')
subplot(1,3,3)
plot(ratios.*theta_SVHJ(vartest),comp3,'r',ratios.*theta_SVHJ(vartest),aux.*min(comp3),'*k')
legend('2nd step-1optm')



%%%Monte-Carlo Routine for SV model
% nmc=100;
% start=zeros(nmc,5);out2=start;out1=start;start2=start;

% for i=1:nmc
% states = sim_states([log(100),0.2^2,5],N,dt,type,theta_SV);%394068); %394068Exemplu de seed foarte bun la 501 instante 776528
% states_estim = [diff(states(:,1)),states(2:end,2)]; %create stationary return series for estimation

% alpha=1e-20;

% theta_SV2=theta_SV*(0.5+rand(1,1));
% start(i,:)=theta_SV2;
% % % % % %     %Find starting values 
% % % % % %     func1 = @(k) critstep1([k,theta_SV2(2:end)],states_estim,dt,type);
% % % % % %     options=optimset('Display','iter');%,'TolX',1e-10,'TolFun',1e-10);
% % % % % %     k=fminsearch(func1,theta_SV2(1),options);
% % % % % %     start2(i,1)=k;
% % % % % %     for j=1:length(theta_SV2)-1
% % % % % %     func1 = @(k) critstep1([theta_SV2(1:j),k,theta_SV2(j+2:end)],states_estim,dt,type);
% % % % % %     options=optimset('Display','iter');%,'TolX',1e-10,'TolFun',1e-10);
% % % % % %     step1=fminsearch(func1,theta_SV2(j+1),options);
% % % % % %     start2(i,j+1)=step1;
% % % % % %     end
% % % % % %     func1 = @(k) critstep1([theta_SV2(1:end-1),k],states_estim,dt,type);
% % % % % %     options=optimset('Display','iter');%,'TolX',1e-10,'TolFun',1e-10);
% % % % % %     k=fminsearch(func1,theta_SV2(end),options);
% % % % % %     start2(i,end)=k;
% % % % % % 
% func1 = @(theta) critstep1(theta,states_estim,dt,type);
% func1(theta_SV)
% options=optimset('Display','iter');%,'TolX',1e-10,'TolFun',1e-10);
% step1=fminsearch(func1,start(i,:),options);
% out1(i,:)=step1;

% [OptMatrix, H0]=funcH0(step1,states_estim, dt, type, alpha);
% func2 = @(theta) critstep2(theta,states_estim,dt,type,H0,OptMatrix);
% func2(theta_SV)
% options=optimset('Display','iter');
% step2=fminsearch(func2,step1,options);
% out2(i,:)=step2;

% [OptMatrix, H0]=funcH0(step2,states_estim, dt, type, alpha);
% func2 = @(theta) critstep2(theta,states_estim,dt,type,H0,OptMatrix);
% func2(theta_SV)
% options=optimset('Display','iter');
% out2(i,:)=fminsearch(func2,step2,options);

% end
% vartest=5;
% comp1=[];
% comp2=[];
% theta_SV3=out2(end,:);
% ratios=0.1:.25:3;
% for rate=ratios
%     theta_SV3(vartest)=theta_SV(vartest).*rate;
%         comp1 = [comp1, critstep1(theta_SV3, states_estim, dt, type)];
%         comp2 = [comp2, critstep2(theta_SV3, states_estim, dt, type, H0, OptMatrix)];
% end
% subplot(1,2,1)
% plot(ratios.*theta_SV(vartest),comp1)
% legend('1st step')
% subplot(1,2,2)
% plot(ratios.*theta_SV(vartest),comp2,'r')
% legend('2nd step')




%  SVHJ Case - Simulate the state vector
% dt=5/250;
% type='SVHJ';
% 
% states = sim_states([log(100),0.2^2,5],N,dt,type,theta_SVHJ,491270);%,871045);%394068); %394068Exemplu de seed foarte bun la 501 instante 776528
% 
% states_estim = [diff(states(:,1)),states(2:end,2:end)]; %create stationary return series for estimation
% 
% alpha=1e-17;
% theta_SVHJ2=theta_SVHJ;
% [OptMatrix, H0]=funcH0(theta_SVHJ2, states_estim, dt, type, alpha);
% out1=[]
% out2=[]
% ratios=0.1:.05:2.1;
% for rate=ratios
% theta_SVHJ=theta_SVHJ2;
% theta_SVHJ(7)=theta_SVHJ2(7).*rate;
% out1=[out1; critstep1(theta_SVHJ, states_estim, dt, type)];
% out2=[out2; critstep2(theta_SVHJ, states_estim, dt, type, H0, OptMatrix)];
% end
% subplot(1,2,1)
% plot(ratios,out1)
% legend('1st step')
% subplot(1,2,2)
% plot(ratios,out2,'r')
% legend('2nd step')

% % % %   
% % % %   tic
% % % %   [deriv_data, true_states] = sim_deriv_data([log(100),0.2^2,5],N,5/250,'SVHJ',theta,[100],[seed]);
% % % %   toc
% % % %   deriv_data(:,7)=deriv_data(:,5);
% % % %   fprintf('\n')
% % % %   %Determine the structure of the data file
% % % %         M=max(deriv_data(:,1));
% % % %         index=1;
% % % %         i=1;    
% % % %         in_day=0;
% % % %         while i<=M
% % % %           TT_marks(i,1)=i;
% % % %           TT_marks(i,2)=0;
% % % %           %Pick the options in each day
% % % %           in_day=1;
% % % %           in_TT=3;
% % % %           TT_now=deriv_data(index,2);
% % % %           while index<=length(deriv_data) && deriv_data(index,1)<=i
% % % %           if(TT_now<deriv_data(index,2))
% % % %               TT_marks(i,in_TT)=in_day-1;
% % % %               in_TT=in_TT+1;
% % % %               TT_now=deriv_data(index,2);
% % % %           end
% % % %           index=index+1;
% % % %           in_day=in_day+1;
% % % %           end  
% % % %           TT_marks(i,in_TT)=in_day-1;
% % % %           i=i+1;
% % % %         end
% % % % 
% % % % %% Estimate parameters based on the simulated sample
% % % %     lb=[0, -0.40, 0.02, 0.5,   2, 0.005, 0.1, -0.8, 0.5, 0.05, 3]; %bound for theta(1) and theta(9) are here bounds for the difference
% % % %     ub=[0.1, -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3, 30, 1.50, 30];
% % % % %     isnotokay=1;
% % % % %     fprintf('Picking theta_start\n');
% % % % %     while isnotokay
% % % % %     thetas=theta+((-1).^round(rand(1,11))).*0.2.*(1-rand(1,11)).*theta;
% % % % %     if sum(imag(boundtofree('SVHJ',thetas,lb,ub)))==0
% % % % %         isnotokay=0;
% % % % %     end
% % % % %     end
% % % % %     fprintf('Theta was picked\n');
% % % % %     theta=thetas;
% % % % %     sim_result=estimate2( deriv_data, TT_marks, 5/250,'SVHJ', theta,lb,ub,theta);
% % % % %     name_for_save=sprintf('deriv_data_%3.0f_%5.0f.mat',N,seed);
% % % % %     save(sprintf('Data Sim/Deriv Samples/%s',name_for_save));
% % % %     
% % % %     
% % % %     
% % % % %     lb=[0, -0.40, 0.02, 0.5,   2, 0.01, 0.1, -0.8,0.1]; %lb SV
% % % % %     ub=[0.07, -0.04, 0.20, 5.0,  10, 0.10, 0.4, -0.3,5 ]; %ub SV
% % % % % theta = lb+rand([1,11]).*(ub-lb);
% % % % % theta(9) = theta(11) + theta(9); %theta(9) stores the difference, theta(11) the level
% % % % % theta(1) = theta(2) + theta(1); %theta(1) stores the difference, theta(2) the level
% % % % test(:,1) = 0.8:0.1:1.5;
% % % % bad_start=theta+0.2.*(1-rand(1,11)).*theta ;   %%%+[0.6.*(1-rand(1,3)).*theta(1:3),zeros([1,8])];
% % % % for q=1:length(test(:,1))
% % % % mod=bad_start;
% % % % mod(9)=theta(9)*test(q,1);
% % % % test(q,2)=quickcrit2( deriv_data, TT_marks, 5/250,'SVHJ', mod,lb,ub,bad_start);
% % % % end
% % % % plot(test(:,1),test(:,2));
% % % % 
% % % % % % % deriv_data(:,4)=deriv_data(:,4)+normrnd(0,0.001,size(deriv_data(:,4)));
% % % % 
% % % % % % 
% % % % % theta1step=theta;
% % % % % X=impl_state_func_fut(deriv_data,TT_marks,'SVHJ',theta);
% % % % % plot([X(:,3),true_states(:,3)])
% % % % %   factor=1e3;
% % % % %   alpha=0.01;
% % % % % H0=funcH0(theta1step, X, 5/250, 'SVHJ', factor);
% % % % % fun = @(theta) estimate_hess(deriv_data, TT_marks, 5/250,'SVHJ',theta',lb,ub,theta1step,X,H0);
% % % % % out=hessian_2sided(fun,theta');

load('chirp');
sound(y);
disp('Done!');  
toc