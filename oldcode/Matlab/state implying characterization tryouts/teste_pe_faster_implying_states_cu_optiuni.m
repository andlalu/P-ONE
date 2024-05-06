clc; clear;
rng(6666);
%% Initialize parameters
kappa = 7; 
vbar = 0.1; 
sigmav = .1; 
rho = -0.5; 
theta_V = [kappa, vbar, sigmav, rho];
eta = 2;
theta_S = eta;
%% State and option price simulation
% Specify the sample details % Small mock-up sample. 
N_timepoints = 2;
true_states = [ 4, 1.5*vbar;
                4, 1.5*vbar;];
% K_vect(1,1,1) = 1;
K_vect = .8:.005:1.2; %moneyness (strike to ATM ratio) 
% K = 1; 
r = 0;
T = 0.01:0.01:1;
theta = [theta_S, theta_V];
[opt_price, opt_price_jacob] = opt_pricing(K_vect,T(1),r,theta,true_states);

%%
% Simulate option prices
tic
% Option pricing parameters:
for i = 1 : length(T)
theta = [theta_S, theta_V];
[opt_price, opt_price_jacob, coeffs] = opt_pricing(K_vect,T(i),r,theta,true_states);
opt_price = reshape(opt_price,[N_timepoints,length(K_vect)]);
for j = 1 : length(K_vect)
 opt_price_table(i,j) = imp_vol_v2_0(exp(true_states(1,1)), T(i), K_vect(j)*exp(true_states(1,1)), opt_price(1,j), r, 'call');
end
end
toc
  cMap = parula(256);
  dataMax = 150;
  dataMin = 1;
  centerPoint = 1;
  scalingIntensity = 4;
  x = 1:length(cMap); 
  x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
  x = scalingIntensity * x/max(abs(x));
  x = sign(x).* exp(abs(x));
  x = x - min(x); x = x*511/max(x)+1; 
  newMap = interp1(x, cMap, 1:512);
  parameter_names = {'\kappa' 'v_{bar}' '\sigma_v' '\rho'};
for parno = 2:length(theta) 
  % Deviate the parameter value
  theta = [theta_S, theta_V];
  theta(parno) = theta(parno)*1.2;
  % Simulate option prices
  tic
  % Option pricing parameters
  for i = 1 : length(T)
   [opt_price, opt_price_jacob, coeffs] = opt_pricing(K_vect,T(i),r,theta,true_states);
   opt_price = reshape(opt_price,[N_timepoints,length(K_vect)]);
   for j = 1 : length(K_vect)
    opt_price_table2(i,j) = imp_vol_v2_0(exp(true_states(1,1)), T(i), K_vect(j)*exp(true_states(1,1)), opt_price(1,j), r, 'call');
   end
   end
  toc  
  % Compute the differences
  differences1 = abs(opt_price_table2-opt_price_table);

  % Deviate the parameter value
  theta = [theta_S, theta_V];
  theta(parno) = theta(parno)*.8;
  % Simulate option prices
  tic
  % Option pricing parameters
  for i = 1 : length(T)
   [opt_price, opt_price_jacob, coeffs] = opt_pricing(K_vect,T(i),r,theta,true_states);
   opt_price = reshape(opt_price,[N_timepoints,length(K_vect)]);
   for j = 1 : length(K_vect)
    opt_price_table2(i,j) = imp_vol_v2_0(exp(true_states(1,1)), T(i), K_vect(j)*exp(true_states(1,1)), opt_price(1,j), r, 'call');
   end
   end
  toc  
  differences = (differences1 + abs(opt_price_table2-opt_price_table))/2;
  
  % Add the contourplot to the figure of all parameters. 
  subplot(floor((length(theta)-1)/2),2,parno-1)
  d = contourf(differences,100,'Linecolor','none'); 
  title(['Implied vol. change due to ' parameter_names(parno-1)]);
  xlabel('Moneyness (K/S in %)');
  ylabel('Maturity (in months)')
  ax = gca;
  ax.XTick = 1:10:length(K_vect);
  ax.XTickLabels = K_vect(1:10:length(K_vect)).*100;
  ax.YTick = [(1:6:length(T)),length(T)];
  ax.YTickLabels = [T(1:6:length(T))*12 T(end)*12];
  c = colorbar;
  c.Label.String = 'Implied vol. change (%)';
  c.Ticks = [0 0.0025 0.005 0.01 0.015 0.02];
  c.TickLabels = {'0%' '0.25%' '0.5%' '1%' '1.5%' '2%'};
  c.Limits = [0 0.02];
  caxis([0 0.02]);
  colormap(newMap)
  hold on;
  [x, y] = find(differences == max(max(differences)));
  plot(x,y,'k+', 'MarkerSize', 14, 'LineWidth', 1.5);
  hold on; 
  txt1 = sprintf('max = %1.3f %s',max(max(differences))*100,'%');
  text(x+2,y+2,txt1,'Interpreter','tex');
  hold off;
end


% theta(2) = 10;
% [~,~,coeffs] = opt_pricing(K_vect,T,r,theta,true_states);
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','CheckGradients',false,'Display','off','SpecifyObjectiveGradient',true);
% minfun = @(state2) fun(state2, K_vect,T,r,theta,true_states(1),opt_price(1,:),coeffs);
% out = lsqnonlin(minfun,1e-2,[],[],options);

% 
% % theta  = [2.46429959524672,0.256909290968590,0.161865908915456,0.0513840767282583,-0.738430244976813];
% tic
% nonISQMLE(theta,N_timepoints,opt_price,K_vect,T,r,states)
% toc
% 
% funcmin_handle = @(theta) -nonISQMLE(theta,N_timepoints,opt_price,K_vect,T,r,states);
% options = optimset('Display','iter','MaxFunEvals',1000);
% %out = fminsearch(funcmin_handle,theta,options);
% out = fmincon(funcmin_handle,theta,[],[],[],[],[1e-5 1e-5, 1e-5, 1e-5, -1],[inf, inf, inf, .5, 1],[],options);
% 
% tic
% nonISQMLE(out,N_timepoints,opt_price,K_vect,T,r,states)
% toc
% %% Setting up a fast state implying routine using the numerical derivatives. 
% 
% function out = nonISQMLE(theta,N_timepoints,opt_price,K_vect,T,r,states)
%   options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','CheckGradients',false,'Display','off','SpecifyObjectiveGradient',true);
%   [~,~,coeffs] = opt_pricing(K_vect,T,r,theta,states);
% % %   implying latent vol from option prices
% %   tic
%   parfor i = 1:N_timepoints
%    minfun = @(state2) fun(state2, K_vect,T,r,theta,states(i,1),opt_price(i,:),coeffs);
%    state2_test(i,2) = max([lsqnonlin(minfun,.2,[],[],options),0]); 
%   end
% %   toc
% % %   computing the QMLE criterion
% %   tic
%   
%   logpdf=zeros(N_timepoints,1);
%   for i=2:N_timepoints
%     [Mu(1),Mu(2),Sigma(1,1),Sigma(2,2),Sigma(1,2)] = heston_qmle_moments(state2_test(i-1,2),theta);
%     Sigma(2,1)=Sigma(1,2); 
%     logpdf(i,1) = logmvnpdf([states(i,1)-states(i-1,1), state2_test(i,2)],Mu,Sigma);
%   end
%   out = sum(logpdf)./(N_timepoints-1);
% %   fprintf('theta: %d %d %d %d %d; crit: %d \n',[theta out]);
%   % toc
% end

% tic
% minfun = @(state2) fun(state2, K_vect,T,r,theta,states(1,1),opt_price(1,:),coeffs);
% state2_test(1,2) = lsqnonlin(minfun,.2,[],[],options);
% for i = 2:N_timepoints
%  minfun = @(state2) fun(state2, K_vect,T,r,theta,states(i,1),opt_price(i,:),coeffs);
%  state2_test(i,2) = lsqnonlin(minfun,state2_test(i-1,2),[],[],options);
% end
% toc

%% QMLE criterion build-up

% % Defining a function that evaluates moments for the QMLE approach. -> Can this be vectorized????
% function [EY, EV, VAR_Y, VAR_V, COV_YV] = heston_qmle_moments(vzero, theta)
% [eta,kappa,vbar,sigmav,rho] = ...
%      deal(theta(:,1),theta(:,2),theta(:,3),theta(:,4),theta(:,5));
% t = 1/252;
% etaSTAR = eta-.5;
% EY = etaSTAR*(t*vbar + ((exp(-kappa*t) - 1)*(vbar - vzero))/kappa);
% EV = vbar - exp(-kappa*t)*(vbar - vzero);
% VAR_V = (sigmav^2*exp(-2*kappa*t)*(exp(kappa*t) - 1)*(2*vzero - vbar + vbar*exp(kappa*t)))/(2*kappa);
% VAR_Y = (2*kappa^2*vzero - 2*kappa^2*vbar + 2*kappa^3*t*vbar + 2*kappa^2*vbar*exp(-kappa*t) - 2*kappa^2*vzero*exp(-kappa*t) + 4*etaSTAR*kappa^2*vzero^2 - 5*etaSTAR^2*sigmav^2*vbar + 2*etaSTAR^2*sigmav^2*vzero - 4*etaSTAR*kappa^2*vzero^2*exp(-kappa*t) + 4*etaSTAR^2*sigmav^2*vbar*exp(-kappa*t) + etaSTAR^2*sigmav^2*vbar*exp(-2*kappa*t) - 2*etaSTAR^2*sigmav^2*vzero*exp(-2*kappa*t) + 2*etaSTAR^2*kappa*sigmav^2*t*vbar - 8*etaSTAR*kappa*rho*sigmav*vbar + 4*etaSTAR*kappa*rho*sigmav*vzero + 8*etaSTAR*kappa*rho*sigmav*vbar*exp(-kappa*t) - 4*etaSTAR*kappa*rho*sigmav*vzero*exp(-kappa*t) + 4*etaSTAR*kappa^2*rho*sigmav*t*vbar + 4*etaSTAR^2*kappa*sigmav^2*t*vbar*exp(-kappa*t) - 4*etaSTAR^2*kappa*sigmav^2*t*vzero*exp(-kappa*t) + 4*etaSTAR*kappa^2*rho*sigmav*t*vbar*exp(-kappa*t) - 4*etaSTAR*kappa^2*rho*sigmav*t*vzero*exp(-kappa*t))/(2*kappa^3);
% COV_YV = (etaSTAR*sigmav^2*vbar + 2*kappa*rho*sigmav*vbar - etaSTAR*sigmav^2*vbar*exp(-2*kappa*t) - 2*etaSTAR*sigmav^2*vzero*exp(-kappa*t) + 2*etaSTAR*sigmav^2*vzero*exp(-2*kappa*t) - 2*kappa*rho*sigmav*vbar*exp(-kappa*t) - 2*etaSTAR*kappa*sigmav^2*t*vbar*exp(-kappa*t) + 2*etaSTAR*kappa*sigmav^2*t*vzero*exp(-kappa*t) - 2*kappa^2*rho*sigmav*t*vbar*exp(-kappa*t) + 2*kappa^2*rho*sigmav*t*vzero*exp(-kappa*t))/(2*kappa^2);
% end


%% fast, efficient log(MVNPDF) function
% function [logp] = logmvnpdf(x,mu,Sigma)
% % outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% % x is NxD, mu is 1xD, Sigma is DxD
% 
% [N,D] = size(x);
% const = -0.5 * D * log(2*pi);
% 
% xc = bsxfun(@minus,x,mu);
% 
% term1 = -0.5 * sum((xc / Sigma) .* xc, 2); % N x 1
% term2 = const - 0.5 * logdet(Sigma);    % scalar
% logp = term1' + term2;
% 
% end
% 
% function y = logdet(A)
% 
% U = chol(A);
% y = 2*sum(log(diag(U)));
% 
% end
% 

%% Option implying function bit. 
function [out1,out2] = fun (state2, K_vect,T,r,theta,state1,opt_price,coeffs)
[out1, out2] = opt_pricing(K_vect,T,r,theta,[state1,state2],coeffs);
out1 = -opt_price + out1;
out2=out2';
end


%%
% % % % % %% Filtering a la christoffersen  - now in a function set-up. 
% % % % % %Parameter set-up for the particle filtering routine
% % % % % theta_filtering = theta;
% % % % % deltat = 1/252;
% % % % % funcmin_handle = @(theta_filtering) -particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call);
% % % % % options = optimset('Display','iter','MaxFunEvals',500);
% % % % % out = fminsearch(funcmin_handle,theta,options);
% % % % % % [criterion_out, v_particles_out] = particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call);
% % % % % %Function:
% % % % % function [criterion,varargout] = particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call)
% % % % %  [eta,kappa,vbar,sigmav,rho] = ...
% % % % %       deal(theta_filtering(:,1),theta_filtering(:,2),theta_filtering(:,3),theta_filtering(:,4),theta_filtering(:,5));
% % % % %  pdf_integral= zeros(N_timepoints-1,1);
% % % % %  ret_particle_PDF = nan(N_timepoints,n_particles);
% % % % %  opt_particle_PDF = nan(N_timepoints,n_particles);
% % % % %  particle_PDF = nan(N_timepoints,n_particles);
% % % % %  coeffs = cf_coeffs(T,theta_filtering,a,b,N_cos);
% % % % %  tic
% % % % %  for k=2:N_timepoints-1
% % % % %    z(1:n_particles/2) = randn(n_particles/2,1);
% % % % %    z(n_particles/2+1:n_particles) = -z(1:n_particles/2);
% % % % %    v_particles(k,:) =  v_particles(k-1,:) + kappa.*(vbar-v_particles(k-1,:))*deltat + sigmav*rho*(states(k,1)-states(k-1,1) - (eta-.5)*v_particles(k-1,:)*deltat) + sigmav*sqrt(v_particles(k-1,:))*sqrt(1-rho^2)*sqrt(deltat).*z;
% % % % %    means = (eta-0.5)*v_particles(k,:)*deltat;
% % % % %    stdevs = sqrt(deltat)*sqrt(v_particles(k,:));
% % % % %    ret_particle_PDF(k,:) = exp(-0.5 * (((states(k+1,1)-states(k,1)) - means)./stdevs).^2)./(sqrt(2*pi).*stdevs);
% % % % %    %Price options for each particle
% % % % %    cf_particles = exp( coeffs(:,1)*ones(1,length(K_vect)*size(v_particles,2)) + (repmat(coeffs(:,2:end)*([states(k,1)*ones(size(v_particles(k,:)))' v_particles(k,:)']'),[1,length(K_vect)]) - 1i*k_grid*pi/(b-a)*(reshape((states(k,1)*ones(length(v_particles(k,:)),1)) + log(K_vect),[1,length(v_particles(k,:))*length(K_vect)]))));
% % % % %    opt_price_particles = exp(-r*T) * exp(reshape((states(k,1)*ones(size(v_particles(k,:)))' + log(K_vect)),[1,length(v_particles(k,:))*length(K_vect)]))...
% % % % %                            .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cf_particles)));
% % % % %    opt_particle_PDF(k,:) = prod(exp(-0.5 * (((reshape(opt_price_particles,length(v_particles(k,:)),length(K_vect))'-opt_price(1,k:N_timepoints:N_timepoints*length(K_vect))')./.01).^2))./(sqrt(2*pi).*.01),1);
% % % % %    particle_PDF(k,:) = ret_particle_PDF(k,:).*opt_particle_PDF(k,:);
% % % % %    pdf_integral(k) = log(mean(particle_PDF(k,:)));
% % % % %    % Re-sample the particle space for using probabilities as weights:
% % % % %    v_particles(k,:) = randsample(v_particles(k,:),n_particles,true,particle_PDF(k,:));
% % % % %  end
% % % % %  toc
% % % % %  criterion = sum(pdf_integral);
% % % % %  disp(theta_filtering);
% % % % %  disp(criterion);
% % % % %  varargout{1} = v_particles; 
% % % % %  varargout{2} = particle_PDF;
% % % % %  varargout{3} = ret_particle_PDF;
% % % % %  varargout{4} = opt_particle_PDF;
% % % % % end

%% Other functions 
% Option pricing function (stand-alone)
function [out,varargout] = opt_pricing(K_vect,T,r,theta,states,varargin)
% Option pricing FFT COS routine coordinates: 
a  = -5.0; %integration lower bound
b  =  5.0; %integration upper bound
c  =    .0; %other integration bound
d  =  b;    %integration bound
N_cos 	= 512;  %Integration bound for COS pricing part
k_grid = cumsum(ones(N_cos,1))-1; %Integration grid values
% Compute the coefficients of the characteristic function
if nargin==5
coeffs = cf_coeffs(T,theta,a,b,N_cos);
else
coeffs = varargin{1};
end
% FFT COS routine set-up:
arg_dma  = k_grid*pi*(d-a)/(b-a);
arg_cma  = k_grid*pi*(c-a)/(b-a);
fact_sin = k_grid*pi/(b-a);
chi = cos(arg_dma)*exp(d) - cos(arg_cma)*exp(c) + fact_sin.*sin(arg_dma)*exp(d) - fact_sin.*sin(arg_cma)*exp(c); 
chi = chi.*((1 + fact_sin.^2).^(-1));
phi = (sin(arg_dma)-sin(arg_cma)).*[0; (b-a).*((k_grid(2:end)*pi).^(-1))] + [(d-c);zeros(N_cos-1,1)];
V_call = 2.0/(b-a)*(chi-phi);
%Sample of option prices:
%Local, fast (vectorized) option pricing routine
cfv2 = exp( coeffs(:,1)*ones(1,length(K_vect)*size(states,1)) + (repmat(coeffs(:,2:end)*states',[1,length(K_vect)]) - 1i*k_grid*pi/(b-a)*(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])))); 
out = exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2)));
%Attempt to compute the first order derivative analytically
if nargout==2
cfv2_fod = exp( coeffs(:,1)*ones(1,length(K_vect)*size(states,1)) + (repmat(coeffs(:,2:end)*states',[1,length(K_vect)]) - 1i*k_grid*pi/(b-a)*(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])))).*coeffs(:,3).*coeffs(:,3); 
varargout{1} = exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2_fod)));
elseif nargout==3
varargout{2} = coeffs;
end
% % % %Approximating the derivative numerically.
% % % h = 1e-5;
% % % states(:,2) = states(:,2) + h;
% % % cfv2_num = exp( coeffs(:,1)*ones(1,length(K_vect)*size(states,1)) + (repmat(coeffs(:,2:end)*states',[1,length(K_vect)]) - 1i*k_grid*pi/(b-a)*(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])))); 
% % % opt_price_num = exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2_num)));
% % % deriv_num = (opt_price_num - opt_price)/h;
% % % states(:,2) = states(:,2) - h;
%'contaminated' option price: 
% opt_price = .1*randn([1,length(states(:,1))*length(K_vect)]) + exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2)));
end

function out = integrate_S_and_V (theta_S, theta_V, logF0, V0, deltat)
  randn1 = randn;
  eta = theta_S(1);
  out(1) = logF0 + (eta-.5)*V0*deltat + sqrt(V0) * randn1*sqrt(deltat);
  k = theta_V(1);
  vbar = theta_V(2);
  sigmav = theta_V(3);
  rho = theta_V(4); %assign parameter values
  out(2) = V0 + k*(vbar-V0)*deltat + sigmav * sqrt(V0) * (rho * randn1 *sqrt(deltat) + sqrt(1-rho^2)* randn* sqrt(deltat));
end
% % % % 
% %% Deprecated 
% % % % function out = integrate_S_and_V_faster (theta_S, theta_V, out_init, deltat, n_steps)
% % % % eta = theta_S(1);
% % % % k = theta_V(1);
% % % % vbar = theta_V(2);
% % % % sigmav = theta_V(3);
% % % % rho = theta_V(4); %assign parameter values
% % % % epsilons = randn(n_steps,2).*sqrt(deltat);
% % % % out(1,:) = out_init;
% % % % for i=2:n_steps+1
% % % %  out(i,2) = out(i-1,2) + k*(vbar-out(i-1,2))*deltat + sigmav * sqrt(out(i-1,2)) * (rho*epsilons(i-1,1) + sqrt(1-rho^2)*epsilons(i-1,2));
% % % %  out(i,1) = out(i-1,1) + (eta-.5).*deltat.*out(i-1,2) + sqrt(out(i-1,2)) * epsilons(i-1,1);
% % % % end
% % % % % out = out(2:n_steps+1,:);
% % % % out = out(n_steps:n_steps+1,:);
% % % % end
% % % % 
% % % % Routine aimed at simulating sv type model under P. 
function out = integrate_V (theta_V, V0, deltat)
  k = theta_V(1);
  vbar = theta_V(2);
  sigmav = theta_V(3);
  rho = theta_V(4); %assign parameter values
  out = V0 + k*(vbar-V0)*deltat + sigmav * sqrt(V0) * (rho * randn *sqrt(deltat) + sqrt(1-rho^2)* randn* sqrt(deltat));
end
