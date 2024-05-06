function out = imply_states_heston (opt_prices, K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho, varargin)
%imply_states Imply latent stochastic states using fully specified AJD model and COS transform based option pricing
% Input: 
%   opt_prices - the panel of option prices
%   K_vect - strikes 
%   T_vect - maturities
%   r_vect - discount rates for each elem in T_vect
%   states_mat - matrix of observable state variable values
%   etaV - affine vol risk premium
%   k - mean reversion 
%   vBar - long run vol state level
%   sigmaV - vol of vol 
%   rho - instantaneous correlation  
% Optional input: 
%   cfCoeffs_mat - matrix of cfCoeffs to be used for pricing. (!)Has to match FFT routine coordinates defined locally. 
% Output: 
%  - (n+1) vector of latent stochastic volatility values
%
% Details: n/a.

%% Checks: no checks implemented.

%% Logic: 

% Optimization settings
options = optimoptions(@lsqnonlin,...
                       'Algorithm','levenberg-marquardt',...
                       'CheckGradients',false,...
                       'Display','off',...                      %'final-detailed',...
                       'SpecifyObjectiveGradient',true,...
                       'OptimalityTolerance', 1.0000e-10,...
                       'FunctionTolerance', 1.0e-10,...
                       'StepTolerance', 1.0e-10);

% Check if ccf coeffs were provided 
if (length(varargin) == 1)
  coeffs = varargin{1};
else 
  [~,~,coeffs] = price_heston_vectorized(K_vect, T_vect, r_vect, NaN, etaV, k, vBar, sigmaV, rho, false, false); 
end
 
% Context variables
N_timepoints = length(states_mat);
N_strikes = length(K_vect);
out = NaN(N_timepoints,1);
% Implying latent vol from option prices
minfun = @(vol) fun(vol, K_vect, T_vect, r_vect, etaV, k, vBar, sigmaV, rho, states_mat(1,1), opt_prices(1:N_timepoints:1+(N_strikes-1)*N_timepoints), coeffs);
out(1) = max([lsqnonlin(minfun, vBar , [], [], options),1e-10]); 
parfor i = 1:N_timepoints
    minfun = @(vol) fun(vol, K_vect, T_vect, r_vect, etaV, k, vBar, sigmaV, rho, states_mat(i,1), opt_prices(i:N_timepoints:i+(N_strikes-1)*N_timepoints), coeffs);
    out(i) = max([lsqnonlin(minfun, vBar, [], [], options),1e-10]); 
end
end

% Auxiliary criterion function for state implying
function [out1, out2] = fun (vol, K_vect, T_vect, r_vect, etaV, k, vBar, sigmaV, rho, states_mat, opt_prices, coeffs)
[out1, out2] = price_heston_vectorized(K_vect, T_vect, r_vect, [states_mat, vol], etaV, k, vBar, sigmaV, rho, false, false, coeffs);
out1 = (-opt_prices(:) + out1(:));%./opt_prices(:);
out2 = out2(:);%./opt_prices(:);
end
