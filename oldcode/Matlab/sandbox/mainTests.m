%Sandbox mode for testing. 
clc; clear;
plots = false;

%checking out the heston state simulation routine:
%seed:
seed = uint8(1);

%parameters: 
etaS = 2; 
k = 7; 
vBar = 0.1; 
sigmaV = .1; 
rho = -0.5; 
etaV = 0;

%coordinates: other
N_sims = 1;
N_timepoints = 252*10;
N_intermediate_draws = 10;
start_LogF0 = log(100);
start_v0 = vBar;
dt = 1/252;

%call: state sim routine
tic 
states_mat = sim_heston(N_sims, N_timepoints, N_intermediate_draws, start_LogF0, start_v0, dt, etaS, etaV, k, vBar, sigmaV, rho, seed);
toc

%plot: first set of states to inspect
if (plots)
  subplot(2,1,1);
  plot(1:1:N_timepoints+1, states_mat(:,1:2:end));
  subplot(2,1,2);
  plot(1:1:N_timepoints+1, states_mat(:,2:2:end));
end

%coordinates: opt pricing coordinates
K_vect = [.9:.05:1.1]; %money-ness levels
T_vect = [.5]; % maturities
r_vect = zeros(size(T_vect)); % discount rates for testing purposes

%parameters: vol risk premium (turn off for testing). 
etaV = 0.0; 

%call: cos transform based option pricing routine vectorized returning call prices
tic 
optPrices_fast = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho);
toc % this bit has to run as fast as possible

%test_out a state implying procedure: 
tic
test = imply_states_heston(optPrices_fast, K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho);
toc

%compare numerical derivative w.r.t. latent vol state vs numerical approximation
function [] = derivative_test
% first test, call out price_heston_vectorized to return derivatives and prices: 
  tic 
  [test0, test2] = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho);
  toc
  
  % second test approximate the derivative numerically 
  states_mat_test = states_mat;
  states_mat_test(:,2) = states_mat_test(:,2) + 0.0001;
  test1 = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat_test, etaV, k, vBar, sigmaV, rho);
  approx_deriv = (test1 - test0)./0.0001;
  deriv_diff = test2 - approx_deriv;
  plot(deriv_diff)
end


%Unit test: Monte Carlo option pricer to check transform prices 
function [] = MC_test()
 MC_test_maturities = [100*dt, 200*dt, 300*dt];
 MC_test_transform_optPrices = price_heston_vectorized(K_vect, MC_test_maturities, zeros(size(MC_test_maturities)), [start_LogF0, start_v0], etaV, k, vBar, sigmaV, rho); 
 % simulation based pricing
 MC_test_NPaths = 30000;
 tolerance = 0.01;  % relative tolerance which if exceeded triggers warning message 
 MC_test_NTimepoints = max(MC_test_maturities./dt);
 MC_test_sim_states_mat = sim_heston(MC_test_NPaths, MC_test_NTimepoints, N_intermediate_draws, start_LogF0, start_v0, dt, 0, etaV, k, vBar, sigmaV, rho, seed);
 MC_test_price = [];
 for tt = MC_test_maturities./dt 
  MC_test_stock_mat = exp(MC_test_sim_states_mat(int32(tt), 1:2:end)); %select log F only
  MC_test_price_row = [];
  for kk = K_vect
   MC_test_price_row = [MC_test_price_row, sum(MC_test_stock_mat(MC_test_stock_mat - kk*exp(start_LogF0)  > 0) - kk*exp(start_LogF0))/MC_test_NPaths];
  end
  MC_test_price = [MC_test_price; MC_test_price_row];
 end
  if any(abs(MC_test_price - MC_test_transform_optPrices)./MC_test_transform_optPrices > tolerance)
   fprintf('Maximum error vs MC too large: %2.3f perc. \n', 100*max(max(abs(MC_test_price - MC_test_transform_optPrices)./MC_test_transform_optPrices)));
  end
end

%Unit test: implying vols from dollar value option prices 
function [] = ImpVol_test()
%call: cos transform based option pricing routine vectorized 
tic 
optPrices = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho);
toc

%call: cos transform based option pricing routine vectorized returning B-S implied vols
tic 
return_implied_vols = true;
use_fast_bs_imp_vol = false;
optVols = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho, return_implied_vols, use_fast_bs_imp_vol);
toc
end