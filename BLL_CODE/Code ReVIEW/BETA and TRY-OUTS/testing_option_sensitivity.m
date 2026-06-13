%% Intro
tic
clc;
clear;
%% Simulation parameters
c=10;
muj0 = -0.05;
mujq0 = -0.15;
sigmaj0 = 0.1;
eta0 = 2.5;
k0 = 5;
vbar0 = 0.01;
sigma0 = 0.20;
rho0 = -0.6;
kl0 = 20; %Starting value is the difference from delta
lambda0 = 0.5;
delta0 = 10;

% % %  SVHJ Case - Price options for varying parameter sets 
dt=5/250;
N=501;
type='SVHJ';
theta_SVHJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ
i=1;
ratios=-.05:.05:.05;
vartest=2;
theta_SVHJ2=theta_SVHJ;
% std_T_graph=(0.1:0.05:1)';
% std_K_graph=(0.8:0.05:1.2)';
% cf_coeff=zeros([4096,4,length(std_T_graph),length(ratios)]);
% vols=zeros([length(std_K_graph),length(std_T_graph)+1,length(ratios)]);
[deriv_data,states,TT_marks]=sim_deriv_data([log(100);0.2^2;5],501,dt,type,theta_SVHJ,1111);
state_imp=zeros([500,3,length(ratios)]);
for rate=ratios
    theta_SVHJ2(vartest)=theta_SVHJ(vartest)+rate
    rate
    state_imp(:,:,i) = impl_state_func_fut(deriv_data,TT_marks,type,theta_SVHJ2);
%     cf_coeff(:,:,:,i)= price_fft_fut(0,0,std_T_graph,0,type,theta_SVHJ2,0,1);
%     vols(:,:,i)= price_fft_fut([log(100);sqrt(0.4);5],std_K_graph.*100 , std_T_graph, 0, type, theta_SVHJ2,1,0,cf_coeff(:,:,:,i));
    i=i+1;
end


%% Plotting the implied state change
figure
set(0,'defaulttextinterpreter','latex')
subplot(2,1,1);
hold all
x=flip(brewermap(length(ratios),'Set1'));
i=1;
tt=find((ratios==0));
for t=1:1:length(ratios)
plot(squeeze(state_imp(:,2,i))-states(2:end,2),'Color',x(i,:),'LineWidth',1);
% if t~=tt
label{i}=['$v$[$\delta$=' num2str((theta_SVHJ(vartest)+ratios(t)),'%2.1f') ']' '$-v^{true}$;~~' ];
% else
% label{i}=['$v^{true}$'];    
% end
i=i+1; 
end
h=title({'Impact of $\delta$ on Implied Latent Volatility',''});
% set(h, 'interpreter', 'latex');
xlabel({'Sample-points'; ''});
ylabel('$v$');
h=legend(label);
set(h, 'interpreter', 'latex');
set(h,'Location','eastoutside','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(2,1,2);
hold all
i=1;
for t=1:1:length(ratios)
plot(squeeze(state_imp(:,3,i))-states(2:end,3),'Color',x(i,:),'LineWidth',1);
% if t~=tt
label{i}=['$\lambda$[$\delta$=' num2str((theta_SVHJ(vartest)+ratios(t)),'%2.1f') ']' '$-\lambda^{true}$;~~' ];
% else
% label{i}=['$v^{true}$'];    
% end
i=i+1; 
end
h=title('Impact of $\delta$ on Implied Latent Jump Intensity');
set(h, 'interpreter', 'latex');
xlabel({'Sample-points'});
ylabel({'$\lambda$'});
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','eastoutside','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

%% Plotting the surface change as a result of parameter change:

% figure
% hold all
% x=flip(brewermap(length(ratios),'RdBu'));
% i=1;
% tt=find((ratios==0));
% for t=1:1:length(ratios)
% if t~=tt
% surf(std_T_graph,std_K_graph,vols(:,2:end,t)-vols(:,2:end,tt),'EdgeColor',x(t,:),'FaceColor',x(t,:));
% set(gcf,'Colormap',brewermap(t,'PuBu'))
% freezeColors
% label{i}=['BSImpVol [$\delta$=' num2str((theta_SVHJ(vartest)+ratios(t)),'%2.1f') ']' ' - BSImpVol [$\delta$=' num2str(theta_SVHJ(vartest),'%2.1f') ']' ];
% i=i+1; 
% end
% 
% if t==tt
% surf(std_T_graph,std_K_graph,vols(:,2:end,t)-vols(:,2:end,tt),'EdgeColor',[.7,.7,.7],'FaceColor',[1,1,1]);
% set(gcf,'Colormap',brewermap(t,'PuBu'))
% freezeColors
% label{i}=['0'];
% i=i+1;
% end
% 
% end
% 
% grid on
% alpha(.8)
% view(-70,15)
% set(gca,'yTick',min(std_K_graph):0.1:max(std_K_graph));
% set(gca,'xTick',min(std_T_graph):0.2:max(std_T_graph))
% set(gca,'zTick',-0.05:0.025:0.05)
% axis([0.1 1 0.8 1.2 -0.05 0.05])
% set(gca,'zTickLabel',{'-5%','-2.5%', '0','2.5%','5%'})
% % set(gca,'yTickLabel',num2str([min(std_K_graph):0.1:max(std_K_graph)].*100,'%2.0f'));
% % set(findall(gcf,'type','axes'),'fontname','Helvetica','fontsize',12) %'fontweight','bold',
% % set(findall(gcf,'type','text'),'fontname','Helvetica','fontsize',12)%'fontweight','bold',
% h=title('Impact of $\delta$ on option prices');
% set(h, 'interpreter', 'latex');
% xlabel({'Maturity (years)'})
% ylabel({'Money-ness (S/K)'})
% zlabel('Implied Volatility');
% h=legend(label);
% set(h, 'interpreter', 'latex')
% set(h,'Location','northwest');
% hold off
% % print -painters -dpdf -r600 surf_k.pdf



% % % Phased out
%         states_imp=impl_state_func_fut(deriv_data, TT_marks, type, theta_SV2);
%         states_estim = [diff(states(:,1)),states_imp(2:end,2)]; %create stationary return series for estimation
%         comp1 = [comp1, critstep1(theta_SV2, states_estim, dt, type)];
%         comp2 = [comp2, critstep2(theta_SV2, states_estim, dt, type, H0, OptMatrix)];
