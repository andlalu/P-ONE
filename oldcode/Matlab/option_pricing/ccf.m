%%% V 3.0.0 <-07.08.2014
function [out] = ccf(u,T,K0,K1,H0,H1,L0,L1,JT)
% Computes the CCF coefficients for AJD process using (a version of) Runge-Kutta

%Call ODE solver
func = @(t,Din) diff(t,Din,K0,K1,H0,H1,L0,L1,JT);
[~,out2] = my_ode45(func,[0;T],[0;1i.*u]);
out(1,:,:)=out2(2:end,:).';
end

%% FUNCTION
function out = diff(t,Din,K0,K1,H0,H1,L0,L1,JT)
Din=Din(2:end);
aux=zeros(length(Din),1);
for i=1:length(Din)
    aux(i,1)=0.5*(Din.')*H1(:,:,i)*Din;
end
out=[K0'*Din + 0.5*(Din.')*H0*Din + L0*JT(Din);K1'*Din  + aux + L1*JT(Din);];
end