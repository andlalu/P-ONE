clc;clear;


% [a, b, c]=meshgrid(-8:2:8,-8:2:8,-8:2:8); x=[a(:),b(:),c(:)];
% w=ones;
% % w=exp(-0.5*sum(x.*x,2));   
% % w=w2./sum(w2);

% [x,w]=nwspgr('GQN', 3, 5);

% x=1:1:1000;
% x=x';
% x=[-x,x,zeros(size(x))];
% [a, b, c]=meshgrid(x,x,x);  x=[a(:),b(:), c(:)];
% x=[zeros(length(x),1),zeros(length(x),1),x'];
% w=exp(-0.5*sum(x.*x,2)/1e6);
% w=w./(sum(w));


% x=[1:1:10];
% x=[x, zeros(size(x)), flip(x);   zeros(size(x)), x, zeros(size(x));     zeros(size(x)),zeros(size(x)),x]';
% w=exp(-0.5*sum(x.*x,2)/1e4);
% w=w./(sum(w));
%% Genz & Keister fully symmetric interpolatory product-rule
[w,x]=fwtpts( 3, 5, 'Norm');
x=x';
w=w';
%%
muj0 = -0.05;
mujq0 = -0.14;
sigmaj0 = 0.06;
eta0 = 2.4;
k0 = 4.8;
vbar0 = 0.01;
sigma0 = 0.22;
rho0 = -0.6;
kl0 = 38; %Starting value is the difference from delta
lambda0 = 0.3;
delta0 = 37.5;

theta = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ

c=1;

muj=theta(1);
mujq=theta(2);
sigmaj=theta(3);
eta=theta(4);
k=theta(5);
vbar=theta(6);
sigma=theta(7);
rho=theta(8);
kl=theta(9);
lambda=theta(10)/c;
delta=theta(11)/c;
%Define CCF system matrixes
K0 = [ 0;    k*vbar;   kl*lambda];
K1 = [ 0,   eta-0.5, -c*(exp(mujq+0.5*sigmaj^2)-1);
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
T=5/250;
% j=1;

for j=1:length(x)
aux=ccf(x(j,:).',T,K0,K1,H00,H1,L0,L1,JT).';
tast(j,1)=exp([1,4,0.015,1/c]*aux);
test(j,1) = exp([4,0.015,1/c]*x(j,:)'.*1i)-exp([1,4,0.015,1/c]*aux);
% test(length(x)-j+1,1)=exp([1,0,0.015,1/c]*conj(aux));
% w(j,1)=exp(-0.5*sum([0,ll,ll]).^2./100);
% j=j+1;
end
% w=w./(sum(w));
subplot(4,2,1)
plot(1:length(x),(tast).*conj(tast))
legend('ccf theta_1')
subplot(4,2,2)
plot(1:length(x),(test).*conj(test).*w)
legend('mom con integrand theta_1')

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

theta = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0]; %theta SVHJ

% c=100;

muj=theta(1);
mujq=theta(2);
sigmaj=theta(3);
eta=theta(4);
k=theta(5);
vbar=theta(6);
sigma=theta(7);
rho=theta(8);
kl=theta(9);
lambda=theta(10)/c;
delta=theta(11)/c;
%Define CCF system matrixes
K0 = [ 0;    k*vbar;   kl*lambda];
K1 = [ 0,   eta-0.5, -c*(exp(mujq+0.5*sigmaj^2)-1);
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
T=5/250;
% j=1;
% range=-10:.5:10;
tic
for j=1:length(x);
aux=ccf(x(j,:).',T,K0,K1,H00,H1,L0,L1,JT).';
tast2(j,1)=exp([1,4,0.015,1/c]*aux);
test2(j,1) = exp([4,0.015,1/c]*x(j,:)'.*1i)-exp([1,4,0.015,1/c]*aux);
% test2(length(x)-j+1,1)=exp([1,0,0.015,1/c]*conj(aux));
% w(j,1)=exp(-0.5*sum([0,ll,ll]).^2./100);
% j=j+1;
end
toc

subplot(4,2,3)
plot(1:length(x),(tast2).*conj(tast2))
legend('ccf theta_2 phi*conj(phi)')
subplot(4,2,4)
plot(1:length(x),(test2).*conj(test2).*w)
legend('mom con integrand theta_2')
subplot(4,2,5)
plot(1:length(x),(tast).*conj(tast)-(tast2).*conj(tast2))
legend('difference between ccf')
subplot(4,2,6)
plot(1:length(x),((tast).*conj(tast)-(tast2).*conj(tast2)).*w)
legend('weighted difference between ccf')
subplot(4,2,7)
plot(1:length(x),(test).*conj(test)-(test2).*conj(test2))
legend('difference between mom cond')
subplot(4,2,8)
plot(1:length(x),((test).*conj(test)-(test2).*conj(test2)).*w)
legend('weighted difference between mom cond')



% % % test=[]
% % % q=-10:.1:10;
% % % q=[zeros(size(q')),zeros(size(q')),q'];
% % % for i=1:length(q)
% % % cf2=ccf(q(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % test(i,1)=(exp(1i*sum(q(i,:).*(X(5,:)),2))-exp(sum(cf2.*[1,0,X(4,2:end)],2)));
% % % end
% % % subplot(3,1,1)
% % % plot(test(:,1).*conj(test(:,1)))
% % % q=-10:.1:10;
% % % q=[zeros(size(q')),q',zeros(size(q'))];
% % % for i=1:length(q)
% % % cf2=ccf(q(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % test(i,2)=(exp(1i*sum(q(i,:).*(X(3,:)),2))-exp(sum(cf2.*[1,0,X(2,2:end)],2)));
% % % end
% % % plot(test(:,2).*conj(test(:,2)))
% % % 
% % % 
% % % q=-1:.1:1;
% % % q=[zeros(size(q')),q',zeros(size(q'))];
% % % for i=1:length(q)
% % % cf2=ccf(q(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % test(i,2)=(exp(1i*sum(q(i,:).*(X(5,:)),2))-exp(sum(cf2.*[1,0,X(4,2:end)],2)));
% % % end
% % % subplot(3,1,2)
% % % plot(test(:,2).*conj(test(:,2)))
% % % 
% % % 
% % % 
% % % q=-10:.1:10;
% % % q=[q',zeros(size(q')),zeros(size(q'))];
% % % for i=1:length(q)
% % % cf2=ccf(q(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % test(i,3)=(exp(1i*sum(q(i,:).*(X(5,:)),2))-exp(sum(cf2.*[1,0,X(4,2:end)],2)));
% % % end
% % % subplot(3,1,3)
% % % plot(test(:,3).*conj(test(:,3)))