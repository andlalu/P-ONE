%% Quick and dirty investigation of self-excitation impact on option prices
std_T_graph=(0.1:0.05:1)';
std_K_graph=(0.8:0.01:1.2)';
%coeff with jump size risk
cf_coeff_save= price_fft_fut(0,0,std_T_graph,0,'SVHJ',theta,0,1);

%coeff without jump size risk 
theta2=theta;
theta2(2)=theta(1);
cf_coeff_save2= price_fft_fut(0,0,std_T_graph,0,'SVHJ',theta2,0,1);




%% Static analysis
ivols = [[0,(0.1:0.05:1)];price_fft_fut([log(100);vbar0;lambda0],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 1, 0,cf_coeff_save)];
ivols2 = [[0,(0.1:0.05:1)];price_fft_fut([log(100);vbar0;lambda0],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta2, 1, 0,cf_coeff_save2)];

prices = price_fft_fut([log(100);vbar0;lambda0],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta, 0, 0,cf_coeff_save);
prices2 = price_fft_fut([log(100);vbar0;lambda0],std_K_graph.*100 , std_T_graph, 0, 'SVHJ', theta2, 0, 0,cf_coeff_save2);

%% Average differences
clearvars d_rel_avg I J d d_rel X Y d_abs_avg
d=prices(:,2:end)-prices2(:,2:end);
d_rel=d./prices(:,2:end);
I=0.925:0.025:1.075;
J=[0.1,0.25,0.5,0.75,1];
for i = 2:length(I)
   for j = 2:length(J)
     [X, Y] =  meshgrid(find((std_K_graph>=I(i-1)).*(std_K_graph<=I(i))==1),find((std_T_graph>=J(j-1)).*(std_T_graph<=J(j))==1));
     d_rel_avg(i-1,j-1)= mean(mean(d_rel([X(:) Y(:)])));
     d_abs_avg(i-1,j-1)= mean(mean(d([X(:) Y(:)])));
   end
end
 d_rel_avg=[I(2:end)', d_rel_avg];
 d_rel_avg=[[0 J(2:end)]; d_rel_avg];
 
 d_abs_avg=[I(2:end)', d_abs_avg];
 d_abs_avg=[[0 J(2:end)]; d_abs_avg];