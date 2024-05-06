function sigma = bs_imp_vol_fwd( obs_price, T, F, k, r)
%BS_IMP_VOL_FWD Implying volatility from call prices using the Black pricing model
%   Inputs:
%    obs_price = observed call price from which volatility should be implied
%    F = forward price for maturity T 
%    K = strike price
%    r = risk free rate used to compute the option price
% 
%   Outputs:
%   out = implied volatility of the call option price

func = @(sigma) bs_price_fwd(exp(-r*T)*obs_price,sigma,T, F, k*F);
    sigma = my_fzero(func,[eps,5]);
end

function out = bs_price_fwd(x,sigma,T, F, K)
%BS_PRICE_DIF Prices undiscounted call option using Black model and computes difference from x
%   Inputs:
%    x = reference value compared to which the difference is determined
%    sigma = volatility level
%    T = option maturity, EXPRESSED IN YEARS
%    S = stock price 
%    K = strike price

% BS d term:
d = (log(F/K) + sigma^2*0.5*T)/(sigma*sqrt(T));

% Note that the normcdf function calculation for a variable x is 0.5 * erfc(-x ./ sqrt(2))
out = x - (0.5 * erfc(-d ./ sqrt(2))*F - 0.5 * erfc(-(d-sigma*sqrt(T)) ./ sqrt(2))*K);

end