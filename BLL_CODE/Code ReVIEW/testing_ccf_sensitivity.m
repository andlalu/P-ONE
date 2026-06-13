%% Intro
clc;clear;

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
lambda0 = 0.5;
delta0 = 10;

theta_SV = [eta0,k0,vbar0,sigma0,rho0]; %theta SV
theta_SVJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,lambda0]; %theta SVJ
theta_SVHJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ

%% Simulate data-set of states & derivatives
dt=5/250;
N=501;
type='SVHJ';
theta_SVHJ = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ
d=1;
ratios=-5:5:5;
vartest=11;
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
    state_imp(:,:,d) = impl_state_func_fut(deriv_data,TT_marks,type,theta_SVHJ2);
% state_imp(:,:,d)=[diff(states),states(2:end,2),states(2:end,3)];
    c=100;
    c2=1;
    theta=theta_SVHJ2;
    theta(10)=theta_SVHJ2(10)./c;
    theta(11)=theta_SVHJ2(11)./c;
     %Allocate parameters    
    muj=theta(1);  
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigma=theta(7);
    rho=theta(8);
    kl=theta(9);
    lambda=theta(10);
    delta=theta(11);
    K0 = [ 0;    k*vbar;   kl*lambda];
    K1 = [ 0,   (eta-0.5), -c*(exp(mujq+0.5*sigmaj^2)-1);
           0,     -k,                           0;
           0,      0,                         -kl; ];
    H00=[ 0,   0,   0;
          0,   0,   0;
          0,   0,   0; ];
    H1(:,:,1) = [ 0,0,0; 0,0,0; 0,0,0;];
    H1(:,:,2) = [ 1,sigma*rho, 0; sigma*rho,sigma^2,0; 0,0,0;];
    H1(:,:,3) = [ 0,0,0; 0,0,0; 0,0,0;];
    L0 = 0;
    L1 = [0;0;c];
    JT = @(beta) (exp(beta(1)*muj+0.5*beta(1).^2*sigmaj^2+delta*beta(3))-1);
    state_imp(:,3,:)=state_imp(:,3,:)./c;
    %% CCF evaluation and plotting routine
    t=-100:.5:100;t=t';
    t1=[t,zeros(size(t)),zeros(size(t))];
    N=size(state_imp,1);
        for i=1:length(t)
        cft1=ccf(t1(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
%         test1a(:,i)=exp(1i*sum(repmat(t1(i,:),[N-1,1]).*(state_imp(2:end,:,d)),2));
%         test1b(:,i)=exp(sum(repmat(cft1,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),state_imp(1:end-1,2:end,d)],2));
        test1(:,i,d)=(exp(1i*sum(repmat(t1(i,:),[N-1,1]).*(state_imp(2:end,:,d)),2))-exp(sum(repmat(cft1,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),state_imp(1:end-1,2:end,d)],2)));
        end
%        
    t=-100:.5:100;t=t';
    t1=[zeros(size(t)),t,zeros(size(t))];
    N=size(state_imp,1);
        for i=1:length(t)
        cft2=ccf(t1(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
        test2(:,i,d)=(exp(1i*sum(repmat(t1(i,:),[N-1,1]).*(state_imp(2:end,:,d)),2))-exp(sum(repmat(cft2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),state_imp(1:end-1,2:end,d)],2)));
        end
% 
    t=-100:.5:100;t=t';
    t1=[zeros(size(t)),zeros(size(t)),t];
    N=size(state_imp,1);
        for i=1:length(t)
        cft3=ccf(t1(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
        test3(:,i,d)=(exp(1i*sum(repmat(t1(i,:),[N-1,1]).*(state_imp(2:end,:,d)),2))-exp(sum(repmat(cft3,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),state_imp(1:end-1,2:end,d)],2)));
        end
d=d+1;  
end 
subplot(3,2,1);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,real(mean(test1(:,:,i))));
label{i }=['$\Re(\overline{MC_1}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_1=\exp(\mathrm{i}[s_1,0,0]\cdot X^\theta) -\phi([s_1,0,0],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_1$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(3,2,2);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,imag(mean(test1(:,:,i))));
label{i}=['$\Im(\overline{MC_1}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_1=\exp(\mathrm{i}[s_1,0,0]\cdot X^\theta) -\phi([s_1,0,0],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_1$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(3,2,3);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,real(mean(test2(:,:,i))));
label{i}=['$\Re(\overline{MC_2}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_2=\exp(\mathrm{i}[0,s_2,0]\cdot X^\theta) -\phi([0,s_2,0],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_2$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(3,2,4);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,imag(mean(test2(:,:,i))));
label{i}=['$\Im(\overline{MC_2}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_2=\exp(\mathrm{i}[0,s_2,0]\cdot X^\theta) -\phi([0,s_2,0],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_2$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(3,2,5);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,real(mean(test3(:,:,i))));
label{i}=['$\Re(\overline{MC_3}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_3=\exp(\mathrm{i}[0,0,s_3]\cdot X^\theta) -\phi([0,0,s_3],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_3$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off

subplot(3,2,6);
hold all
i=1;
for g=1:1:length(ratios)
plot(t,imag(mean(test3(:,:,i))));
label{i}=['$\Im(\overline{MC_3}[\delta=$' num2str((theta_SVHJ(vartest)+ratios(g)),'%2.2f') '$])$' ];
i=i+1; 
end
h=title('Sample Average of $MC_3=\exp(\mathrm{i}[0,0,s_3]\cdot X^\theta) -\phi([0,0,s_3],X^\theta,\Delta;\theta)$','FontSize',9);
set(h, 'interpreter', 'latex');
h=xlabel('$s_3$');
set(h, 'interpreter', 'latex');
h=legend(label);
set(h, 'interpreter', 'latex')
set(h,'Location','northeast','fontname','Helvetica','orientation','vertical');
set(findall(gcf,'type','axes'),'fontname','Helvetica');
hold off



