function out = cf_coeffs(T, theta,a,b,N_cos)  
%%% V0.1 (alpha) 23-01-2017. European option pricing routine using the COS approach by Fang & Osterlee (2008) 
%%% A. Lalu (a.lalu@uva.nl)
%% AJD system matrixes (a la DPS2000) for CCF computations:
 % This case was especially designed for the 
%  len=3;
 %Allocate parameters
 [eta,k,vbar,sigmav,rho] = ...
     deal(theta(:,1),theta(:,2),theta(:,3),theta(:,4),theta(:,5));
 % Ad-hoc, Heston specific analytic solution for the coefficients of the characteristic function
 u = ((cumsum(ones(N_cos,1))-1)*pi/(b-a))*1i; %integration grid coordinates
 bb = u.*sigmav*rho-k;
 aa = u.*(1-u);
 gamma = sqrt(bb.^2+aa*sigmav.^2);
 out(:,1) = -k*vbar*((gamma+bb)/sigmav^2*T + 2/sigmav^2*log(1-(gamma+bb)./(2*gamma) .* (1-exp(-gamma*T))));
 out(:,2) = u;
 out(:,3) = - aa.*(1-exp(-gamma*T))./(2*gamma - (gamma+bb).*(1-exp(-gamma*T)));
    
% % %     
% % % %  %Define CCF system matrixes
% % %  K0  = [ 0; k.*vbar;];
% % %  K1  = [ 0,-0.5;
% % %          0,-k; ];
% % %  H00 = zeros(2,2);
% % %  H1(:,:,1) = zeros(2,2);
% % %  H1(:,:,2) = [ 1, sigmav*rho; sigmav*rho, sigmav^2;];
% % %  L0 = 0;
% % %  L1 = [0; 0;];
% % %  JT = @(beta) (0);
% % %  %Allocate function handle
% % %  char_func = @(arg) ccf([arg;0;],T,K0,K1,H00,H1,L0,L1,JT);
% % % 
% % % %% Option pricing routine:
% % % k = cumsum(ones(N_cos,1))-1; %Integration grid values
% % % %Function returns ccf coefficients
% % % out2=nan([N_cos,len,length(T)]);
% % %  parfor i=1:N_cos
% % %          out2(i,:,:) = char_func(k(i)*pi/(b-a));
% % %  end
        
end

