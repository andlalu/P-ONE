% %% Quick and dirty investigation of self-excitation impact on option prices
clear;
clc;
load('data_withstates.mat');
std_T_graph=[0.1 0.5 1]';
std_K_graph=[0.95 1 1.05]';
%coeff with jump size risk
cf_coeff_save= price_fft_fut(0,0,std_T_graph,0,'SVHJ',theta,0,1);
lambda2=(kl0*lambda0)/(kl0-delta0);
theta2=[muj0 mujq0 sigmaj0 eta0 k0 vbar0 sigma0 rho0 lambda2];
cf_coeff_save2= price_fft_fut(0,0,std_T_graph,0,'SVJ',theta2,0,1);

%% Lambda Path generation
T=0.4; %e.g. 1 year impact analysis?????????
delta=5/250; %say weekly observations of the intensity process
dt=(5/250)/100;
nt=floor(T/delta);
% Simulate the intensity process path to know what we are dealing with
U=[1;zeros(floor(T/dt)-1,1)];
lambda=[lambda0+delta0;zeros(floor(T/dt)-1,1)];
for i=2:floor(T/dt)
lambda(i) = lambda(i-1)+ kl0*(lambda0-lambda(i-1))*dt;% + delta0*(U(i-1)<lambda(i-1)*dt);
end
        lambdalong=lambda;
        lambda=lambda(1:100:end);
        nsims=10;
        d=zeros(length(std_K_graph),length(std_T_graph),nt);
        d_rel=zeros(length(std_K_graph),length(std_T_graph),nt);
for sim=1:1000
        %% Vol Path generation
        v=[vbar0;zeros(floor(T/dt)-1,1)];
        W = sqrt(dt).*randn([T/dt,2]);
        s=[log(100);zeros(floor(T/dt)-1,1)];
        s2=[log(100);zeros(floor(T/dt)-1,1)];
        for i=2:T/dt 
        v(i) = max(v(i-1)+k0*(vbar0-v(i-1))*dt + sigma0*sqrt(v(i-1))*(rho0*W(i-1,1)+sqrt(1-rho0^2)*W(i-1,2)),0);
        s(i) = s(i-1) +((-exp(muj0+0.5*sigmaj0^2)+1)*lambdalong(i-1)+(eta0-0.5)*v(i-1))*dt + sqrt(v(i-1))*W(i-1,1) + -0.1*U(i-1);
        s2(i) = s2(i-1)+((-exp(muj0+0.5*sigmaj0^2)+1)*lambda2   +(eta0-0.5)*v(i-1))*dt + sqrt(v(i-1))*W(i-1,1) + -0.1*U(i-1);
        end
        v=v(1:100:end);
        s=s(1:100:end);
        s2=s2(1:100:end);
        % 


        %% Dynamic analysis
        for t=1:nt

        % vols(:,:,t)= price_fft_fut([s(t);v(t);lambda(t)],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 1, 0,cf_coeff_save);
        % vols2(:,:,t)= price_fft_fut([s2(t);v(t)],std_K_graph.*100 , std_T_graph, 0, 'SVJ', theta2, 1, 0,cf_coeff_save2);
        prices(:,:,t) = price_fft_fut([s(t);v(t);lambda(t)],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 0, 0,cf_coeff_save);
        prices2(:,:,t) = price_fft_fut([s(t);v(t)],std_K_graph.*100 , std_T_graph, 0, 'SVJ', theta2, 0, 0,cf_coeff_save2);

        % clearvars I J  X Y d_rel d
        d(:,:,t)= d(:,:,t)+prices(:,2:end,t)-prices2(:,2:end,t);
        d_rel(:,:,t)=d_rel(:,:,t)+(prices(:,2:end,t)-prices2(:,2:end,t))./prices(:,2:end,t);
        % I=0.925:0.025:1.075;
        % J=[0.1,0.25,0.5,0.75,1];
        % for i = 2:length(I)
        %    for j = 2:length(J)
        %      [X, Y] =  meshgrid(find((std_K_graph>=I(i-1)).*(std_K_graph<=I(i))==1),find((std_T_graph>=J(j-1)).*(std_T_graph<=J(j))==1));
        %      d_rel_avg(i-1,j-1,t)= mean(mean(d_rel([X(:) Y(:)])));
        %      d_abs_avg(i-1,j-1,t)= mean(mean(d([X(:) Y(:)])));
%         end
% % %  d_rel_avg=[I(2:end)', d_rel_avg];
% % %  d_rel_avg=[[0 J(2:end)]; d_rel_avg];
% % %  
% % %  d_abs_avg=[I(2:end)', d_abs_avg];
% % %  d_abs_avg=[[0 J(2:end)]; d_abs_avg];
        end
end 
d=d./sim;
d_rel=d_rel./sim;

% % % % % % % % % Average differences
% % % % % % % % % clearvars d_rel_avg I J d d_rel X Y d_abs_avg
% % % % % % % % % d=prices(:,2:end,1)-prices2(:,2:end,1);
% % % % % % % % % d_rel=d./prices(:,2:end,1);
% % % % % % % % % I=0.925:0.025:1.075;
% % % % % % % % % J=[0.1,0.25,0.5,0.75,1];
% % % % % % % % % for i = 2:length(I)
% % % % % % % % %    for j = 2:length(J)
% % % % % % % % %      [X, Y] =  meshgrid(find((std_K_graph>=I(i-1)).*(std_K_graph<=I(i))==1),find((std_T_graph>=J(j-1)).*(std_T_graph<=J(j))==1));
% % % % % % % % %      d_rel_avg(i-1,j-1)= mean(mean(d_rel([X(:) Y(:)])));
% % % % % % % % %      d_abs_avg(i-1,j-1)= mean(mean(d([X(:) Y(:)])));
% % % % % % % % %    end
% % % % % % % % % end
% % % % % % % % %  d_rel_avg=[I(2:end)', d_rel_avg];
% % % % % % % % %  d_rel_avg=[[0 J(2:end)]; d_rel_avg];
% % % % % % % % %  
% % % % % % % % %  d_abs_avg=[I(2:end)', d_abs_avg];
% % % % % % % % %  d_abs_avg=[[0 J(2:end)]; d_abs_avg];5