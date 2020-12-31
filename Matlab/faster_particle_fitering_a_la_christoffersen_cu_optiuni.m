clc; clear;
rng(6667);
%% Initialize parameters
kappa = 7; 
vbar = 0.1; 
sigmav = .1; 
rho = -0.5; 
theta_V = [kappa, vbar, sigmav, rho];
eta = 2;
theta_S = eta;
%% State and option price simulation
% Specify the sample details
N_timepoints = 1000;
true_states(:,1) = [ 4, vbar];
% K_vect(1,1,1) = 1;
K_vect = [1]; %moneyness (strike to ATM ratio) 
K = 1; 
r = 0;
T = 1/12;
% Simulate state series
tic
for j= 2:N_timepoints*100
true_states(:,j) = integrate_S_and_V (theta_S, theta_V, true_states(1,j-1), true_states(2,j-1), 1/(252*100));
end
states = true_states(:,1:100:N_timepoints*100)';
toc
% Simulate option prices
tic
% Option pricing parameters:
theta = [theta_S, theta_V];
% Option pricing FFT COS routine coordinates: 
a  = -5.0; %integration lower bound
b  =  5.0; %integration upper bound
c  =  .0; %other integration bound
d  =  b;    %integration bound
N_cos 	= 512;  %Integration bound for COS pricing part
k_grid = cumsum(ones(N_cos,1))-1; %Integration grid values
% Compute the coefficients of the characteristic function
coeffs = cf_coeffs(T,theta,a,b,N_cos);
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
opt_price = exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2)));
%'contaminated' option price: 
% opt_price = .1*randn([1,length(states(:,1))*length(K_vect)]) + exp(-r*T) * exp(reshape((states(:,1) + log(K_vect)),[1,length(states(:,1))*length(K_vect)])) .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cfv2)));
toc

%% Particle initialisation
tic 
% Trebui sa initializam un nr de particole sa vedem ce iese.
 %  - Assume particles are initialized by choosing stopped process simulated values. 
% for each particle 
n_particles = 100;
v_old = ones(1,n_particles)*vbar;
v_new = v_old;%nan(1,n_particles);
parfor j = 1 : n_particles 
  for i=1:100
    v_new(j) = integrate_V(theta_V, v_old(j), 1/(252));
    v_old(j) = v_new(j); 
  end
end
v_particles(1,:) = v_new;
toc
%% Filtering a la christoffersen  - now in a function set-up. 
%Parameter set-up for the particle filtering routine
theta_filtering = theta;
deltat = 1/252;
funcmin_handle = @(theta_filtering) -particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call);
options = optimset('Display','iter','MaxFunEvals',1000);
out = fminsearch(funcmin_handle,theta,options);
% [criterion_out, v_particles_out] = particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call);
%Function:
function [criterion,varargout] = particle_filtering(theta_filtering,N_timepoints,n_particles,v_particles,deltat,states,opt_price,K_vect,r,T,a,b,N_cos,k_grid,V_call)
 [eta,kappa,vbar,sigmav,rho] = ...
      deal(theta_filtering(:,1),theta_filtering(:,2),theta_filtering(:,3),theta_filtering(:,4),theta_filtering(:,5));
 pdf_integral= zeros(N_timepoints-1,1);
 ret_particle_PDF = nan(N_timepoints,n_particles);
 opt_particle_PDF = nan(N_timepoints,n_particles);
 particle_PDF = nan(N_timepoints,n_particles);
 coeffs = cf_coeffs(T,theta_filtering,a,b,N_cos);
 tic
 for k=2:N_timepoints-1
   z(1:n_particles/2) = randn(n_particles/2,1);
   z(n_particles/2+1:n_particles) = -z(1:n_particles/2);
   v_particles(k,:) =  v_particles(k-1,:) + kappa.*(vbar-v_particles(k-1,:))*deltat + sigmav*rho*(states(k,1)-states(k-1,1) - (eta-.5)*v_particles(k-1,:)*deltat) + sigmav*sqrt(v_particles(k-1,:))*sqrt(1-rho^2)*sqrt(deltat).*z;
   means = (eta-0.5)*v_particles(k,:)*deltat;
   stdevs = sqrt(deltat)*sqrt(v_particles(k,:));
   ret_particle_PDF(k,:) = exp(-0.5 * (((states(k+1,1)-states(k,1)) - means)./stdevs).^2)./(sqrt(2*pi).*stdevs);
   %Price options for each particle
   cf_particles = exp( coeffs(:,1)*ones(1,length(K_vect)*size(v_particles,2)) + (repmat(coeffs(:,2:end)*([states(k,1)*ones(size(v_particles(k,:)))' v_particles(k,:)']'),[1,length(K_vect)]) - 1i*k_grid*pi/(b-a)*(reshape((states(k,1)*ones(length(v_particles(k,:)),1)) + log(K_vect),[1,length(v_particles(k,:))*length(K_vect)]))));
   opt_price_particles = exp(-r*T) * exp(reshape((states(k,1)*ones(size(v_particles(k,:)))' + log(K_vect)),[1,length(v_particles(k,:))*length(K_vect)]))...
                           .* real(sum(bsxfun(@times,[.5;ones(N_cos-1,1)].*exp(-1i*(k_grid*a*pi/(b-a))).*V_call,cf_particles)));
   opt_particle_PDF(k,:) = prod(exp(-0.5 * (((reshape(opt_price_particles,length(v_particles(k,:)),length(K_vect))'-opt_price(1,k:N_timepoints:N_timepoints*length(K_vect))')./.01).^2))./(sqrt(2*pi).*.01),1);
   particle_PDF(k,:) = ret_particle_PDF(k,:).*opt_particle_PDF(k,:);
   pdf_integral(k) = log(mean(particle_PDF(k,:)));
   % Re-sample the particle space for using probabilities as weights:
   v_particles(k,:) = randsample(v_particles(k,:),n_particles,true,particle_PDF(k,:));
 end
 toc
 criterion = sum(pdf_integral);
 disp(theta_filtering);
 disp(criterion);
 varargout{1} = v_particles; 
 varargout{2} = particle_PDF;
 varargout{3} = ret_particle_PDF;
 varargout{4} = opt_particle_PDF;
end

%% Other functions 
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
