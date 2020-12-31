function [ states_mat ] = sim_heston(nSims, n, m, logF0, v0, deltaT, etaS, etaV, k, vBar, sigmaV, rho, varargin)
%sim_heston Simulates Heston DGP state vector sample paths. 
% Input: 
%   nSims - # of state sample paths to simulate
%   n - # of timepoints to sim per path
%   m - # of increments making up each timepoint
%   logF0 - initial log forward state 
%   v0 - initial local vol state
%   deltaT - time increment between points
%   etaS - affine diffusive risk premium
%   etaV - affine vol risk premium
%   k - mean reversion 
%   vBar - long run vol state level
%   sigmaV - vol of vol 
%   rho - instantaneous correlation  
% Optional input: 
%   seed - number used for the seed of the random number generator 
% Output: 
%  - (n+1) x 2 matrix of state vector values
%
% Details: simulates the Heston DGP state values following Euler discretization of the SDE system: 
% $$d\log(F_t) = (eta_S - 1/2) V_t dt + \sqrt{V_t} dW_{t,1}$$
% $$dV_t = k(V_{bar} - V_t)dt + eta_V V_t dt + \sigma_V \sqrt{V_t} (\rho dW_{t,1} + \sqrt{1- \rho^2} dW_{t,2}$$

%% Checks: no checks implemented.

%% Logic: draws brownian increments and loops to integrate over them. 

%Varargin seed
  if ~isempty(varargin) && isinteger(varargin{1}) && varargin{1} > 0
   rng(varargin{1});
  end

%Draw brownian increments
  dW = randn(n*m, 2*nSims) .* sqrt(deltaT/m);
  F_index = 1:2:2*nSims-1;
  V_index = 2:2:2*nSims;
  dW(:, V_index) = rho * dW(:, F_index) + sqrt(1-rho*rho) * dW(:, V_index);
  
%Compute the state vector values 
  states_mat = NaN(n+1, 2*nSims);
  states_mat(1, :) = repmat([logF0, v0], [1,nSims]);
  bufferStates = NaN(m+1, 2*nSims);
  for i = 1:n
   bufferStates(1, :) = states_mat(i, :);
   for j = 2:m+1
    bufferStates(j, F_index) = bufferStates(j-1, F_index) + (etaS-.5) * deltaT/m * bufferStates(j-1, V_index) + ...
                                + sqrt(bufferStates(j-1, V_index)) .* dW((i-1) * m + (j-1), F_index);
    bufferStates(j, V_index) = bufferStates(j-1, V_index) * (1 + etaV*deltaT/m) + k * (vBar - bufferStates(j-1, V_index)) * deltaT/m + ... 
                                + sigmaV * sqrt(bufferStates(j-1, V_index)) .* dW((i-1) * m + j-1, V_index);
   end
   states_mat(i+1, :) = bufferStates(end, :);
  end
  
%Housekeeping
  clearvars -except states_mat
  
end

