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
N_timepoints = 252*5;
tic
true_states(:,1) = [ 4, vbar];
for j= 2:N_timepoints*100
true_states(:,j) = integrate_S_and_V (theta_S, theta_V, true_states(1,j-1), true_states(2,j-1), 1/(252*100));
end
toc
states = true_states(:,1:100:N_timepoints*100)';
% states(1,:) = [0  diff(states(1,:))];
%% Particle initialisation
tic 
% Trebui sa initializam un nr de particole sa vedem ce iese.
 %  - Assume particles are initialized by choosing stopped process simulated values. 
% for each particle 
n_particles = 500;
v_old = ones(1,n_particles)*vbar;
v_new = v_old;%nan(1,n_particles);
% parfor j = 1 : n_particles 
%   for i=1:100
%     v_new(j) = integrate_V(theta_V, v_old(j), 1/(252));
%     v_old(j) = v_new(j); 
%   end
% end
v_particles(1,:) = v_new;
 toc
%% Filtering 
% sa incercam filtering pentru prima perioada 
% propagam fiecare particula
no_sims = 10; %number of intermediate Euler simulation scheme particle simulation
sim_state = nan(N_timepoints,2,n_particles);
pdf_integral= zeros(N_timepoints-1,1);
particle_PDF = nan(N_timepoints,n_particles);
tic
for k=2:N_timepoints
 s_start = states(k-1,1);
 s_fin = states(k,1);
 subv_particles = v_particles(k-1,:);
 for j=1:n_particles 
  for i=1:10
    aux = integrate_S_and_V_faster (theta_S, theta_V, [s_start,subv_particles(j)], 1/(252*(no_sims)),no_sims);
    td_mean = (eta-0.5)*aux(1,2)*1/(252*(no_sims)) + (aux(2,2)-aux(1,2) - kappa*(vbar-aux(1,2))*1/(252*(no_sims)))/(rho*sigmav);
    td_std = sqrt(1/(252*(no_sims)))*sqrt(1-rho^2)*sqrt(aux(1,2))/sqrt(rho^2);
    aux2(i) = exp(-0.5 * (((s_fin-aux(1,1)) - td_mean)./td_std).^2) ./ (sqrt(2*pi) .* td_std);
  end  
    particle_PDF(k,j) = mean(aux2(i));
 v_particles(k,j) = aux(2,2);
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
% Routine aimed at simulating sv type model under P. 
function out = integrate_V (theta_V, V0, deltat)
  k = theta_V(1);
  vbar = theta_V(2);
  sigmav = theta_V(3);
  rho = theta_V(4); %assign parameter values
  out = V0 + k*(vbar-V0)*deltat + sigmav * sqrt(V0) * (rho * randn *sqrt(deltat) + sqrt(1-rho^2)* randn* sqrt(deltat));
end

function out = integrate_S_and_V_faster (theta_S, theta_V, out_init, deltat, n_steps)
eta = theta_S(1);
k = theta_V(1);
vbar = theta_V(2);
sigmav = theta_V(3);
rho = theta_V(4); %assign parameter values
epsilons = randn(n_steps,2).*sqrt(deltat);
out(1,:) = out_init;
for i=2:n_steps+1
 out(i,2) = out(i-1,2) + k*(vbar-out(i-1,2))*deltat + sigmav * sqrt(out(i-1,2)) * (rho*epsilons(i-1,1) + sqrt(1-rho^2)*epsilons(i-1,2));
 out(i,1) = out(i-1,1) + (eta-.5).*deltat.*out(i-1,2) + sqrt(out(i-1,2)) * epsilons(i-1,1);
end
% out = out(2:n_steps+1,:);
out = out(n_steps:n_steps+1,:);
end


%% Old code
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
