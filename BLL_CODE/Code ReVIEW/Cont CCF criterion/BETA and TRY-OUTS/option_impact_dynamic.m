%% Quick and dirty investigation of self-excitation impact on option prices
std_T_graph=(0.1:0.05:1)';
std_K_graph=(0.8:0.01:1.2)';
%coeff with jump size risk
cf_coeff_save= price_fft_fut(0,0,std_T_graph,0,'SVHJ',theta,0,1);

%% Path generation
T=1; %e.g. 1 year impact analysis?????????
delta=5/250; %say weekly observations of the intensity process
dt=(5/250)/100;
nt=floor(T/delta);
% Simulate the intensity process path to know what we are dealing with
U=[1;zeros(floor(T/dt)-1,1)];
lambda=[lambda0+delta0;zeros(floor(T/dt)-1,1)];
for i=2:floor(T/dt)
lambda(i) = lambda(i-1)+ kl0*(lambda0-lambda(i-1))*dt;% + delta0*(U(i-1)<lambda(i-1)*dt);
end
lambda=lambda(1:100:end);

%% Dynamic analysis
for t=1:nt
vols(:,:,t)= price_fft_fut([log(100);vbar0;lambda(t)],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 1, 0,cf_coeff_save);
prices(:,:,t) = price_fft_fut([log(100);vbar0;lambda(t)],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 0, 0,cf_coeff_save);
end
% %% Average differences
% clearvars d_rel_avg I J d d_rel X Y d_abs_avg
% d=prices(:,2:end)-prices2(:,2:end);
% d_rel=d./prices(:,2:end);
% I=0.925:0.025:1.075;
% J=[0.1,0.25,0.5,0.75,1];
% for i = 2:length(I)
%    for j = 2:length(J)
%      [X, Y] =  meshgrid(find((std_K_graph>=I(i-1)).*(std_K_graph<=I(i))==1),find((std_T_graph>=J(j-1)).*(std_T_graph<=J(j))==1));
%      d_rel_avg(i-1,j-1)= mean(mean(d_rel([X(:) Y(:)])));
%      d_abs_avg(i-1,j-1)= mean(mean(d([X(:) Y(:)])));
%    end
% end
%  d_rel_avg=[I(2:end)', d_rel_avg];
%  d_rel_avg=[[0 J(2:end)]; d_rel_avg];
%  
%  d_abs_avg=[I(2:end)', d_abs_avg];
%  d_abs_avg=[[0 J(2:end)]; d_abs_avg];