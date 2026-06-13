function sigma = bs_imp_vol_fut( obs_price, T, F, K, r)
%BS_IMP_VOL Implying volatility from call prices using the Black pricing model for options written on futures
%   Inputs:
%    obs_price = observed call price from which volatility should be implied
%    F = Futures price for maturity T 
%    K = strike price
%    r = risk free rate used to compute the option price
%    q = dividend rate (set to 0 if not considered)
% 
%   Outputs:
%   out = implied volatility of the call option price

func = @(sigma) bs_price_fut(obs_price,sigma,T, F, K, r);
    sigma = my_fzero(func,[eps,3]);
end

function out = bs_price_fut(x,sigma,T, F, K, r)
%BS_PRICE_DIF Prices call option using Black-Scholes model and computes difference from x
%   Inputs:
%    x = reference value compared to which the difference is determined
%    sigma = volatility level
%    T = option maturity, EXPRESSED IN YEARS
%    S = stock price 
%    K = strike price
%    r = risk free rate used to compute the option price
%    q = dividend rate 
d = (log(F/K) + sigma^2*0.5*T)/(sigma*sqrt(T));


% Note that the normcdf function calculation for a variable x is 0.5 * erfc(-x ./ sqrt(2))

out = x - exp(-r*T)*(0.5 * erfc(-d ./ sqrt(2))*F - 0.5 * erfc(-(d-sigma*sqrt(T)) ./ sqrt(2))*K);

end