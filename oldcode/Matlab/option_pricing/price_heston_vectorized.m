function [out, varargout] = price_heston_vectorized(K_vect, T_vect, r_vect, states_mat, etaV, k, vBar, sigmaV, rho, varargin)
%price_heston_vectorized Efficiently prices European options using COS transform a la Fang & Osterlee (2008) 
% Input: 
%   K_vect - strikes 
%   T_vect - maturities
%   r_vect - discount rates for each elem in T_vect
%   states_mat - matrix of state variable values
%   etaV - affine vol risk premium
%   k - mean reversion 
%   vBar - long run vol state level
%   sigmaV - vol of vol 
%   rho - instantaneous correlation  
% Optional input: 
%   cfCoeffs_mat - matrix of cfCoeffs to be used for pricing. (!)Has to match FFT routine coordinates defined locally. 
% Output: 
%  - (n+1) x 2 matrix of state vector values
%
% Details: n/a.

%% Checks: no checks implemented.

%% Logic: 

% Option pricing FFT COS routine coordinates: 
 a  = -5.0;      %integration lower bound
 b  =  5.0;      %integration upper bound
 c  =   .0;      %other integration bound
 d  =    b;      %integration bound
 N_cos 	= 256;   %integration bound for COS pricing part
 k_grid = cumsum(ones(N_cos,1))-1; %Integration grid values
 lenK = length(K_vect);
 lenS = size(states_mat, 1);
 lenT = length(T_vect);
 widS = size(states_mat, 2);
 
% Compute the coefficients of the characteristic function
switch length(varargin)
 case 0
   u = ((cumsum(ones(N_cos,1))-1)*pi/(b-a))*1i; %integration grid coordinates
   coeffs = cf_coeffs_heston(u, T_vect, etaV, k, vBar, sigmaV, rho);
   return_implied_vols = false;
   use_fast_bs_imp_vol = false;
 case 1
   u = ((cumsum(ones(N_cos,1))-1)*pi/(b-a))*1i; %integration grid coordinates
   coeffs = cf_coeffs_heston(u, T_vect, etaV, k, vBar, sigmaV, rho);
   return_implied_vols = varargin{1};
   use_fast_bs_imp_vol = false;
 case 2
  u = ((cumsum(ones(N_cos,1))-1)*pi/(b-a))*1i; %integration grid coordinates
  coeffs = cf_coeffs_heston(u, T_vect, etaV, k, vBar, sigmaV, rho);
  return_implied_vols = varargin{1}; 
  use_fast_bs_imp_vol = varargin{2};
 case 3
   return_implied_vols = varargin{1}; 
   use_fast_bs_imp_vol = varargin{2};
   coeffs = varargin{3};
end

%quick return of coeffs
if isnan(states_mat)
 out = [];
 varargout{1} = [];
 varargout{2} = coeffs;
 return 
end


% FFT COS routine set-up:
 arg_dma  = k_grid*pi*(d-a)/(b-a);
 arg_cma  = k_grid*pi*(c-a)/(b-a);
 fact_sin = k_grid*pi/(b-a);
 chi = cos(arg_dma)*exp(d) - cos(arg_cma)*exp(c) + fact_sin.*sin(arg_dma)*exp(d) - fact_sin.*sin(arg_cma)*exp(c); 
 chi = chi.*((1 + fact_sin.^2).^(-1));
 phi = (sin(arg_dma)-sin(arg_cma)).*[0; (b-a).*((k_grid(2:end)*pi).^(-1))] + [(d-c);zeros(N_cos-1,1)];
 V_call = 2.0/(b-a)*(chi-phi);

%TODO: Adjust for states_mat containing multiple samples
%e.g.: 
 if widS > 2 
  error('price_heston_vectorized: Multiple samples not supported.'); 
 end

%Sample of option prices using fast routine (vectorized over strike levels K_vect):
 aux_logK = log(K_vect);
 aux_logKreshaped = reshape(aux_logK(ones(1, lenS),:), [], 1)'; 
 aux_coeffCorr  = - 1i * k_grid * pi / (b-a);
 aux_sumWeights = [.5; ones(N_cos-1,1)];
 aux_expansionRow = ones(1, lenK * lenS);
 out = NaN(lenT, lenK * lenS);
 out2 = NaN(lenT, lenK*lenS);

 for tt=1:lenT
  
  % compute the cf for cos transform arg values
  cf = exp(repmat(coeffs(:, 1, tt) + coeffs(:, 3:end, tt) * states_mat(:, 2:end)',  [1, lenK]) ... % (betas * latent states) repeated for each strike
           + aux_coeffCorr * (aux_logKreshaped  + a)); % cos correction term for each N_cos and strike 
  
  % compute the option price
  out(tt,:) = exp(-r_vect(tt) * T_vect(tt) + aux_logKreshaped + reshape(repmat(states_mat(:, 1),[1, lenK]), [], 1)')...
                .* real(sum(((V_call .* aux_sumWeights) *  aux_expansionRow) .* cf, 1));
               
  % compute the gradient w.r.t. latent state if required (!) only works for Heston
  if nargout==2
  out2(tt,:) = exp(-r_vect(tt) * T_vect(tt) + aux_logKreshaped + reshape(repmat(states_mat(:, 1),[1, lenK]), [], 1)')...
                .* real(sum(((V_call .* aux_sumWeights) *  aux_expansionRow) .* cf .* coeffs(:, 3:end, tt), 1));
  end
  
  % return b-s implied vols if required
  if return_implied_vols && use_fast_bs_imp_vol
   for col_i = 1:lenK*lenS
      out(tt,col_i) = fast_bs_imp_vol_fwd(out(tt,col_i), T_vect(tt), exp(states_mat(rem(col_i-1, lenS)+1,1)), K_vect(fix((col_i-1)/lenS)+1), r_vect(tt));
   end
  elseif return_implied_vols
    for col_i = 1:lenK*lenS
      out(tt,col_i) = bs_imp_vol_fwd(out(tt,col_i), T_vect(tt), exp(states_mat(rem(col_i-1, lenS)+1,1)), K_vect(fix((col_i-1)/lenS)+1), r_vect(tt));
    end
  end
 end
 
 % assign additional return variables
 if nargout==2 
   varargout{1} = out2;
 elseif nargout==3
   varargout{2} = coeffs;
 end
 
 function [out_tensor] = cf_coeffs_heston(u, T_vect, etaV, k, vBar, sigmaV, rho)
  %cfCoeffs_heston Returns a tensor containing the CCF coefficients for a Heston model needed for COS pricing. 
  
   % Ad-hoc, Heston specific analytic solution for the coefficients of the characteristic function
   vBar_star = k/(k - etaV) * vBar;
   k_star = k - etaV;
   bb = u.*sigmaV*rho-k_star;
   aa = u.*(1-u);
   gamma = sqrt(bb.^2+aa*sigmaV.^2);
   out_tensor(:,1,:) = -k_star*vBar_star*((gamma+bb)/sigmaV^2*T_vect + 2/sigmaV^2*log(1-(gamma+bb)./(2*gamma) .* (1-exp(-gamma*T_vect))));
   out_tensor(:,2,:) = u * ones(1, length(T_vect));  %repmat(u,[1,length(T_vect)]);
   out_tensor(:,3,:) = - aa.*(1-exp(-gamma*T_vect))./(2*gamma - (gamma+bb).*(1-exp(-gamma*T_vect)));
   
   % Housekeeping local function
   clearvars vBar_star k_star bb aa gamma
  end
 
%Housekeeping: 
 clearvars -except out varargout
 
end