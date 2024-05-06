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
%% State simulation
% Specify the sample details
N_timepoints = 252;
true_states(:,1) = [ 4, vbar];
% Simulate state series
tic
for j= 2:N_timepoints*100
true_states(:,j) = integrate_S_and_V (theta_S, theta_V, true_states(1,j-1), true_states(2,j-1), 1/(252*100));
end
states = true_states(:,1:100:N_timepoints*100)';
% Simulate option prices
theta = [theta_S, theta_V];
coeffs = price_cos_fut([0;0], 0, 3/12, 0, 'SV', theta, 0, 1);
parfor j = 1 : N_timepoints 
 opt_price(j,:) = price_cos_fut(states(j,:)', exp(states(j,1)), 3/12, 0, 'SV', theta, 0, 0, coeffs);
end
toc
%% Particle initialisation
tic 
% Trebui sa initializam un nr de particole sa vedem ce iese.
 %  - Assume particles are initialized by choosing stopped process simulated values. 
% for each particle 
n_particles = 500;
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
%% Filtering a la christoffersen 
% sa incercam filtering pentru prima perioada 
% propagam fiecare particula
pdf_integral= zeros(N_timepoints-1,1);
particle_PDF = nan(N_timepoints,n_particles);
deltat = 1/252;
tic
for k=2:N_timepoints-1
  for j=1:2:n_particles    
    % Primul pas sa recuperam inovatia din progresia lui S_t la S_t+1.
    z = randn;
    v_particles(k,j) = v_particles(k-1,j) + kappa*(vbar-v_particles(k-1,j))*deltat + sigmav*rho*(states(k,1)-states(k-1,1) - (eta-.5)*v_particles(k-1,j)*deltat) + sigmav*sqrt(v_particles(k-1,j))*sqrt(1-rho^2)*sqrt(deltat)*z;
    v_particles(k,j+1) = v_particles(k-1,j+1) + kappa*(vbar-v_particles(k-1,j+1))*deltat + sigmav*rho*(states(k,1)-states(k-1,1) - (eta-.5)*v_particles(k-1,j+1)*deltat) - sigmav*sqrt(v_particles(k-1,j+1))*sqrt(1-rho^2)*sqrt(deltat)*z;
    % After propagating compute the mean and variance of the pdf
    td_mean1 = (eta-0.5)*v_particles(k,j)*deltat;
    td_std1 = sqrt(deltat)*sqrt(v_particles(k,j));
    opt_price1(1,:) = price_cos_fut([states(k,1);v_particles(k,j)], exp(states(k,1)), 3/12, 0, 'SV', theta, 0, 0, coeffs);
    particle_PDF(k,j) = exp(-0.5 * (((states(k+1,1)-states(k,1)) - td_mean1)./td_std1).^2) ./ (sqrt(2*pi) .* td_std1) * exp(-0.5 * ((opt_price(k,2)-opt_price1(1,2))./.001).^2)./(sqrt(2*pi).*.001);
    td_mean2 = (eta-0.5)*v_particles(k,j+1)*deltat;
    td_std2 = sqrt(deltat)*sqrt(v_particles(k,j+1));
    opt_price2(1,:) = price_cos_fut([states(k,1);v_particles(k,j+1)], exp(states(k,1)), 3/12, 0, 'SV', theta, 0, 0, coeffs);
    particle_PDF(k,j+1) = exp(-0.5 * (((states(k+1,1)-states(k,1)) - td_mean2)./td_std2).^2) ./ (sqrt(2*pi) .* td_std2) * exp(-0.5 * ((opt_price(k,2)-opt_price2(1,2))./.001).^2)./(sqrt(2*pi).*.001);
  end
 pdf_integral(k) = log(mean(particle_PDF(k,:)));
 % Re-sample the particle space for using probabilities as weights:
 v_particles(k,:) = randsample(v_particles(k,:),n_particles,true,particle_PDF(k,:));
end
toc
criterion = sum(pdf_integral);
disp(criterion);
% Speed-up routine:vectorize draws + antithetic draws. 
 
%% Functions 
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

%% Deprecated 
% % % function out = integrate_S_and_V_faster (theta_S, theta_V, out_init, deltat, n_steps)
% % % eta = theta_S(1);
% % % k = theta_V(1);
% % % vbar = theta_V(2);
% % % sigmav = theta_V(3);
% % % rho = theta_V(4); %assign parameter values
% % % epsilons = randn(n_steps,2).*sqrt(deltat);
% % % out(1,:) = out_init;
% % % for i=2:n_steps+1
% % %  out(i,2) = out(i-1,2) + k*(vbar-out(i-1,2))*deltat + sigmav * sqrt(out(i-1,2)) * (rho*epsilons(i-1,1) + sqrt(1-rho^2)*epsilons(i-1,2));
% % %  out(i,1) = out(i-1,1) + (eta-.5).*deltat.*out(i-1,2) + sqrt(out(i-1,2)) * epsilons(i-1,1);
% % % end
% % % % out = out(2:n_steps+1,:);
% % % out = out(n_steps:n_steps+1,:);
% % % end

% % % % Routine aimed at simulating sv type model under P. 
function out = integrate_V (theta_V, V0, deltat)
  k = theta_V(1);
  vbar = theta_V(2);
  sigmav = theta_V(3);
  rho = theta_V(4); %assign parameter values
  out = V0 + k*(vbar-V0)*deltat + sigmav * sqrt(V0) * (rho * randn *sqrt(deltat) + sqrt(1-rho^2)* randn* sqrt(deltat));
end
