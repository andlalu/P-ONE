function sigma = imp_vol_v2_0(underlying, tau, strike, price, rate, type)


%%% ORIGINAL DISCLAIMER & DETAILS 
% About LetsBeRational:
% ~~~~~~~~~~~~~~~~~~~~~
% 
% The source code of LetsBeRational resides at www.jaeckel.org/LetsBeRational.7z .
% 
% ======================================================================================
% Copyright © 2013-2014 Peter Jaeckel.
% 
% Permission to use, copy, modify, and distribute this software is freely granted,
% provided that this notice is preserved.
% 
% WARRANTY DISCLAIMER
% The Software is provided "as is" without warranty of any kind, either express or implied,
% including without limitation any implied warranties of condition, uninterrupted use,
% merchantability, fitness for a particular purpose, or non-infringement.
% ======================================================================================

% from __future__ import division, print_function

% from math import fabs, sqrt, log, exp

% from py_lets_be_rational.numba_helper import maybe_jit
% from py_lets_be_rational.exceptions import BelowIntrinsicException, AboveMaximumException

% from py_lets_be_rational.constants import *
% from py_lets_be_rational.rationalcubic import *

% from py_lets_be_rational.erf_cody import erfcx_cody

% from py_lets_be_rational.normaldistribution import inverse_norm_cdf
% from py_lets_be_rational.normaldistribution import norm_cdf
% from py_lets_be_rational.normaldistribution import norm_pdf

DBL_EPSILON = eps;
SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
FOURTH_ROOT_DBL_EPSILON = sqrt(SQRT_DBL_EPSILON);
EIGHTH_ROOT_DBL_EPSILON = sqrt(FOURTH_ROOT_DBL_EPSILON);
SIXTEENTH_ROOT_DBL_EPSILON = sqrt(EIGHTH_ROOT_DBL_EPSILON);
DBL_MIN = realmin;
DBL_MAX = realmax;
SQRT_DBL_MIN = sqrt(DBL_MIN);
SQRT_DBL_MAX = sqrt(DBL_MAX);

VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC = -DBL_MAX;
VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM = DBL_MAX;

ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
SQRT_TWO_PI = 2.506628274631000502415765284811045253006986740610;

TWO_PI = 6.283185307179586476925286766559005768394338798750;
SQRT_PI_OVER_TWO = 1.253314137315500251207882642405522626503493370305;  % sqrt(pi/2) to avoid misinterpretation.
SQRT_THREE = 1.732050807568877293527446341505872366942805253810;
SQRT_ONE_OVER_THREE = 0.577350269189625764509148780501957455647601751270;
TWO_PI_OVER_SQRT_TWENTY_SEVEN = 1.209199576156145233729385505094770488189377498728;  % 2*pi/sqrt(27)
PI_OVER_SIX = 0.523598775598298873077107230546583814032861566563;

implied_volatility_maximum_iterations = 3;
asymptotic_expansion_accuracy_threshold = -10;
small_t_expansion_of_normalized_black_threshold = 2 * SIXTEENTH_ROOT_DBL_EPSILON;

if strcmp(type,'call')
  q=1;
elseif strcmp(type,'put')
  q=-1;
else 
  warning('Unrecognised option type');
end
sigma = nan(length(price),1);
for i = 1:length(price)
sigma(i,1) = implied_volatility_from_a_transformed_rational_guess(price(i), underlying, strike(i), tau, q, implied_volatility_maximum_iterations);
end

end

function out = householder_factor(newton, halley, hh3)
out = (1 + 0.5 * halley * newton) / (1 + newton * (halley + hh3 * newton / 6));
end

function [f,fp,fpp] = compute_f_lower_map_and_first_two_derivatives(x, s)
TWO_PI = 6.283185307179586476925286766559005768394338798750;
SQRT_THREE = 1.732050807568877293527446341505872366942805253810;
SQRT_ONE_OVER_THREE = 0.577350269189625764509148780501957455647601751270;
TWO_PI_OVER_SQRT_TWENTY_SEVEN = 1.209199576156145233729385505094770488189377498728;  % 2*pi/sqrt(27)
PI_OVER_SIX = 0.523598775598298873077107230546583814032861566563;

    ax = abs(x);
    z = SQRT_ONE_OVER_THREE * ax / s;
    y = z * z;
    s2 = s * s;
    Phi = normcdf(-z);
    phi = normpdf(z);
    fpp = PI_OVER_SIX * y / (s2 * s) * Phi * (8 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8) - 8 * x * x) * Phi / phi) * exp(2 * y + 0.25 * s2);
    if is_below_horizon(s)
        fp = 1;
        f = 0;
    else
        Phi2 = Phi * Phi;
        fp = TWO_PI * y * Phi2 * exp(y + 0.125 * s * s);
        if is_below_horizon(x)
            f = 0;
        else
            f = TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (Phi2 * Phi);
        end
    end
end

function  [f,fp,fpp] = compute_f_upper_map_and_first_two_derivatives(x, s)
SQRT_PI_OVER_TWO = 1.253314137315500251207882642405522626503493370305; 
    f = normcdf(-0.5 * s);
    if is_below_horizon(x)
        fp = -0.5;
        fpp = 0;
    else
        w = square(x / s);
        fp = -0.5 * exp(0.5 * w);
        fpp = SQRT_PI_OVER_TWO * exp(w + 0.125 * s * s) * w / s;
    end
end

function out = square(x)
    out =  x*x;
end

function out = inverse_f_lower_map(x, f)
 SQRT_THREE = 1.732050807568877293527446341505872366942805253810;
 TWO_PI_OVER_SQRT_TWENTY_SEVEN = 1.209199576156145233729385505094770488189377498728;
 if is_below_horizon(f)
   out = 0;
 else
   out = abs(x / (SQRT_THREE * norminv(power(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * abs(x)), 1. / 3.))));
 end
end

function out = inverse_f_upper_map(f)
  out =  -2. * norminv(f);
end

function out = is_below_horizon(x)
    % Set this to 0 if you want positive results for (positive) denormalized inputs, else to DBL_MIN.
    % Note that you cannot achieve full machine accuracy from denormalized inputs!
    DENORMALIZATION_CUTOFF = 0;
  %This weeds out denormalized (a.k.a. 'subnormal') numbers.
  out =  abs(x) < DENORMALIZATION_CUTOFF ;
end

  function out = normalized_black_call_using_norm_cdf(x, s)
    h = x / s;
    t = 0.5 * s;
    b_max = exp(0.5 * x);
    b = normcdf(h + t) * b_max - normcdf(h - t) / b_max;
    out =  abs(max(b, 0.0));
  end


function out = asymptotic_expansion_of_normalized_black_call(h, t)
ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
    e=(t/h)*(t/h);
    r=((h+t)*(h-t));
    q=(h/r)*(h/r);
    intermediate = -6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))));
    asymptotic_expansion_sum = 2.0+q*intermediate;
    b = ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*(t/r)*asymptotic_expansion_sum;
    out = abs(max(b , 0.));
end

function out = small_t_expansion_of_normalized_black_call(h, t)
ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
SQRT_TWO_PI = 2.506628274631000502415765284811045253006986740610;
    a = 1+h*(0.5*SQRT_TWO_PI)*erfcx(-ONE_OVER_SQRT_TWO*h);
    w=t*t;
    h2=h*h;
    expansion = 2*t*(a+w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*(-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*(-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))));
    b = ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*expansion;
    out =  abs(max(b,0.0));
end

function out  = normalised_black_call_using_erfcx(h, t)
    ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
    b = 0.5 * exp(-0.5*(h*h+t*t)) * ( erfcx(-ONE_OVER_SQRT_TWO*(h+t)) - erfcx(-ONE_OVER_SQRT_TWO*(h-t)) );
    out =  abs(max(b,0.0));
end

function s = unchecked_normalised_implied_volatility_limited_iterations(beta, x, q,  N)
    DBL_EPSILON = eps;
    % Subtract intrinsic.
    if q*x > 0
        beta = abs(max(beta-normalised_intrinsic(x, q), 0.));
        q = -q;
    end
    % Map puts to calls
    if q < 0
        x = -x;
        q = -q;
    end
    if beta <= 0 % For negative or zero prices we return 0.
        s = 0;
    end
    % Set this to 0 if you want positive results for (positive) denormalized inputs, else to DBL_MIN.
    % Note that you cannot achieve full machine accuracy from denormalized inputs!
    DENORMALIZATION_CUTOFF = 0;
    if beta < DENORMALIZATION_CUTOFF % For positive but denormalized (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
       s = 0;
    end
    b_max = exp(0.5*x);
    if beta >= b_max
        warning('AboveMaximumException');
    end
    iterations = 0;
    direction_reversal_count = 0 ;
    DBL_MAX = realmax;
    DBL_MIN = realmin;
    f = -DBL_MAX;
    s = -DBL_MAX; 
    ds = s; 
    ds_previous = 0;
    s_left = DBL_MIN;
    s_right = DBL_MAX;
    % The temptation is great to use the optimised form b_c = exp(x/2)/2-exp(-x/2)Â·Phi(sqrt(-2Â·x)) but that would require implementing all of the above types of round-off and over/underflow handling for this expression, too.
    s_c = sqrt(abs(2*x));
    b_c = normalised_black_call(x, s_c);
    v_c = normalised_vega(x, s_c);
    % Four branches.
    if beta < b_c
        s_l = s_c - b_c/v_c;
        b_l = normalised_black_call(x,s_l);
        if beta < b_l
            [f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2] = compute_f_lower_map_and_first_two_derivatives(x, s_l);
            r_ll = convex_rational_cubic_control_parameter_right_side(0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,d2_f_lower_map_l_d_beta2,true);
            f = rational_cubic_interpolation(beta,0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,r_ll);
            if ~(f > 0)  % This can happen due to roundoff truncation for extreme values such as |x|>500.
                % We switch to quadratic interpolation using f(0)â?¡0, f(b_l), and f'(0)â?¡1 to specify the quadratic.
                t = beta/b_l;
                f = (f_lower_map_l*t + b_l*(1-t)) * t;
            end
            s = inverse_f_lower_map(x, f);
            s_right = s_l;
            DBL_EPSILON = eps;
            while and(iterations < N, abs(ds) > DBL_EPSILON * s)
                if ds*ds_previous < 0
                    direction_reversal_count = direction_reversal_count + 1;
                end
                if and(iterations>0,any([3==direction_reversal_count,any([s<s_left,s>s_right])])) %Believe this is the correct interpretation of the logical 
                   % If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                   % NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                   s = 0.5*(s_left+s_right);
                   if (s_right-s_left) <= DBL_EPSILON*s
                        break
                   direction_reversal_count = 0;
                   ds = 0;
                   end
                end
                ds_previous = ds;
                b = normalised_black_call(x,s);
                bp = normalised_vega(x, s);
                if and(b>beta,s<s_right)
                    s_right = s;
                elseif and(b<beta,s>s_left)
                    s_left = s; % Tighten the bracket if applicable.
                end
                if or(b<=0,bp<=0) % Numerical underflow. Switch to binary nesting for this iteration.
                   ds = 0.5*(s_left+s_right)-s;
                else
                   ln_b=log(b);
                   ln_beta=log(beta);
                   bpob=bp/b;
                   h=x/s;
                   b_halley = h*h/s-s/4;
                   newton = (ln_beta-ln_b)*ln_b/ln_beta/bpob;
                   halley = b_halley-bpob*(1+2/ln_b);
                   b_hh3 = b_halley*b_halley- 3 * square(h / s) - 0.25;
                   hh3 = b_hh3+ 2 * square(bpob) * (1 + 3 / ln_b * (1 + 1 / ln_b)) - 3 * b_halley * bpob * (1 + 2 / ln_b);
                   ds = newton * householder_factor(newton, halley, hh3);
                end
                ds = max(-0.5*s , ds);
                s = s + ds;
                iterations = iterations + 1;
%                 return  %% HERE??
            end
        else
            v_l = normalised_vega(x, s_l);
            r_lm = convex_rational_cubic_control_parameter_right_side(b_l,b_c,s_l,s_c,1/v_l,1/v_c,0.0,false);
            s = rational_cubic_interpolation(beta,b_l,b_c,s_l,s_c,1/v_l,1/v_c,r_lm);
            s_left = s_l;
            s_right = s_c;
        end
    else
        if v_c > DBL_MIN
          s_h = s_c+(b_max-b_c)/v_c;
        end
        b_h = normalised_black_call(x,s_h);
        if beta <= b_h
            v_h = normalised_vega(x, s_h);
            r_hm = convex_rational_cubic_control_parameter_left_side(b_c,b_h,s_c,s_h,1/v_c,1/v_h,0.0,false);
            s = rational_cubic_interpolation(beta,b_c,b_h,s_c,s_h,1/v_c,1/v_h,r_hm);
            s_left = s_c;
            s_right = s_h;
        else
            DBL_MAX = realmax;
            SQRT_DBL_MAX = sqrt(DBL_MAX);
            [f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2] = compute_f_upper_map_and_first_two_derivatives(x, s_h);
            if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX < SQRT_DBL_MAX
                r_hh = convex_rational_cubic_control_parameter_left_side(b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,d2_f_upper_map_h_d_beta2,true);
                f = rational_cubic_interpolation(beta,b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,r_hh);
            end
            if f <= 0
                h=b_max-b_h;
                t=(beta-b_h)/h;
                f = (f_upper_map_h*(1-t) + 0.5*h*t) * (1-t); % We switch to quadratic interpolation using f(b_h), f(b_max)â?¡0, and f'(b_max)â?¡-1/2 to specify the quadratic.
            end
            s = inverse_f_upper_map(f);
            s_left = s_h ;
            if beta > 0.5*b_max % Else we better drop through and let the objective function be g(s) = b(x,s)-beta.                
                while and(iterations < N,abs(ds) > DBL_EPSILON * s)
                    if ds*ds_previous < 0
                        direction_reversal_count = direction_reversal_count + 1;
                    end
                    if and(iterations>0,or( 3==direction_reversal_count,~and(s>s_left,s<s_right)))
                        % If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                        % NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                        s = 0.5*(s_left+s_right);
                    end
                    if (s_right-s_left<=DBL_EPSILON*s)
                        break
                    end
                    direction_reversal_count = 0;
                    ds = 0;
                    ds_previous=ds;
                    b = normalised_black_call(x,s);
                    bp = normalised_vega(x, s);
                    if and(b>beta,s<s_right)
                        s_right=s;
                    elseif and(b<beta, s>s_left)
                        s_left=s; %Tighten the bracket if applicable.
                    end
                    if or(b>=b_max,bp<=DBL_MIN) % Numerical underflow. Switch to binary nesting for this iteration.
                        ds = 0.5*(s_left+s_right)-s;
                    else
                        b_max_minus_b = b_max-b;
                        g = log((b_max-beta)/b_max_minus_b);
                        gp = bp/b_max_minus_b;
                        b_halley = square(x / s) / s - s / 4;
                        b_hh3 = b_halley*b_halley- 3 * square(x / (s * s)) - 0.25;
                        newton = -g/gp;
                        halley = b_halley+gp;
                        hh3 = b_hh3+gp*(2*gp+3*b_halley);
                        ds = newton * householder_factor(newton, halley, hh3);
                    ds = max(-0.5*s, ds);
                    s = s + ds;
                    iterations = iterations + 1;
                    return %%HERE??
                    end
                end
            end
        end
    end 
    while and(iterations < N,abs(ds) > DBL_EPSILON * s)
        if ds*ds_previous < 0
            direction_reversal_count = direction_reversal_count + 1;
        end
        if and(iterations>0,or(3==direction_reversal_count,~and(s>s_left,s<s_right)))
            % If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
            % NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
            s = 0.5*(s_left+s_right);
            if (s_right-s_left) <= DBL_EPSILON*s
                break
            end
            direction_reversal_count = 0;
            ds = 0;
        end
        ds_previous = ds;
        b = normalised_black_call(x, s);
        bp = normalised_vega(x, s);
        if and(b > beta,s < s_right)
            s_right = s;
        elseif and(b < beta, s > s_left)
            s_left = s; %Tighten the bracket if applicable.
        newton = (beta-b)/bp;
        halley = square(x / s) / s - s / 4;
        hh3 = halley*halley- 3 * square(x / (s * s)) - 0.25;
        ds = max(-0.5 * s, newton * householder_factor(newton, halley, hh3));
        s = s + ds;
        iterations = iterations + 1;
        end
        return
    end
end
    
function out =  normalised_vol_transf_rational_guess_lim(beta, x, q, N)
    
    % Map in-the-money to out-of-the-money
    if q * x > 0
        beta = beta - normalised_intrinsic(x, q);
        q = -q;
    end
    if beta < 0
       warning('BelowIntrinsicException')
    end
    out = unchecked_normalised_implied_volatility_limited_iterations(beta, x, q, N);
end

function out = implied_volatility_from_a_transformed_rational_guess_with_lim(price, F, K, T, q, N)
    if q < 0 
      intrinsic = abs(max(K-F,0));
    else 
      intrinsic = abs(max(F-K,0));
    end
    if price < intrinsic
        warning('BelowIntrinsicException');
    end
    if q < 0
      max_price = K;
    else 
      max_price = F;
    end
    if price >= max_price
        warning('AboveMaximumException');
    end
    x = log(F / K);
    % Map in-the-money to out-of-the-money
    if q * x > 0
        price = abs(max(price - intrinsic, 0.0));
        q = -q;
    end
    out =  unchecked_normalised_implied_volatility_limited_iterations(price / (sqrt(F) * sqrt(K)), x, q, N) / sqrt(T);
end


function out = nor_imp_vol_from_a_transf_rational_guess(beta, x, q) 
    out =  normalised_vol_transf_rational_guess_lim(beta, x, q, implied_volatility_maximum_iterations);
end

function out = implied_volatility_from_a_transformed_rational_guess(price, F, K, T, q, implied_volatility_maximum_iterations)
    out =  implied_volatility_from_a_transformed_rational_guess_with_lim(price, F, K, T, q, implied_volatility_maximum_iterations);
end

function out = normalised_vega(x, s)
DBL_MIN = realmin;
SQRT_DBL_MIN = sqrt(DBL_MIN);
ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
    ax = abs(x);
    if ax <= 0
        out = ONE_OVER_SQRT_TWO_PI * exp(-0.125 * s * s);
    else
       if or(s <= 0,s <= ax * SQRT_DBL_MIN) 
         out = 0;
       else
         out = ONE_OVER_SQRT_TWO_PI * exp(-0.5 * (square(x / s) + square(0.5 * s)));
       end
    end
end

function out = normalised_intrinsic(x, q)
    if q * x <= 0
      out = 0;
      return
    end
    x2 = x * x;
    DBL_EPSILON = eps;
SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
FOURTH_ROOT_DBL_EPSILON = sqrt(SQRT_DBL_EPSILON);
    if x2 < 98 * FOURTH_ROOT_DBL_EPSILON % The factor 98 is computed from last coefficient: â^sâ^s92897280 = 98.1749
        out =  abs(max((sign(q))*x * (1 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))), 0.0 ) );
        return
    end
    b_max = exp(0.5 * x);
    one_over_b_max = 1 / b_max;
    out = abs(max(sign(q)*(b_max - one_over_b_max), 0.));
end

function out = normalised_intrinsic_call(x)
   out = normalised_intrinsic(x, 1);
end

function  out = normalised_black_call(x, s)
    if x > 0
        out = normalised_intrinsic_call(x) + normalised_black_call(-x, s);
        return
    end
    ax = abs(x);
    % Set this to 0 if you want positive results for (positive) denormalized inputs, else to DBL_MIN.
    % Note that you cannot achieve full machine accuracy from denormalized inputs!
    DENORMALIZATION_CUTOFF = 0;
    if s <= ax * DENORMALIZATION_CUTOFF
        out = normalised_intrinsic_call(x);
        return
    end
    % Denote h := x/s and t := s/2. We evaluate the condition |h|>|Î·|, i.e., h<Î·  &&  t < Ï"+|h|-|Î·|  avoiding any
    % divisions by s , where Î· = asymptotic_expansion_accuracy_threshold  and Ï" =
    % small_t_expansion_of_normalized_black_threshold .
    DBL_EPSILON = eps;
    SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
    FOURTH_ROOT_DBL_EPSILON = sqrt(SQRT_DBL_EPSILON);
    EIGHTH_ROOT_DBL_EPSILON = sqrt(FOURTH_ROOT_DBL_EPSILON);
    SIXTEENTH_ROOT_DBL_EPSILON = sqrt(EIGHTH_ROOT_DBL_EPSILON);
    asymptotic_expansion_accuracy_threshold = -10;
    small_t_expansion_of_normalized_black_threshold = 2 * SIXTEENTH_ROOT_DBL_EPSILON;
    if and(x < s * asymptotic_expansion_accuracy_threshold, 0.5 * s * s + x < s * (small_t_expansion_of_normalized_black_threshold + asymptotic_expansion_accuracy_threshold))
        % Region 1.
        out = asymptotic_expansion_of_normalized_black_call(x / s, 0.5 * s);
        return
    end
    if 0.5 * s < small_t_expansion_of_normalized_black_threshold
        % Region 2.
        out = small_t_expansion_of_normalized_black_call(x / s, 0.5 * s);
        return
    end
    % When b is more than, say, about 85% of b_max=exp(x/2), then b is dominated by the first of the two terms in the
    %  Black formula, and we retain more accuracy by not attempting to combine the two terms in any way. We evaluate
    % the condition h+t>0.85  avoiding any divisions by s.
    if x + 0.5 * s * s > s * 0.85
        % Region 3.
        out = normalized_black_call_using_norm_cdf(x, s);
        return
    end
    % Region 4.
    out = normalised_black_call_using_erfcx(x / s, 0.5 * s);
end

function out = normalised_black(x, s, q)
  out =  normalised_black_call(sign(q).*x, s);  % Reciprocal-strike call-put equivalence
end


function out =  black(F, K, sigma, T, q)
    intrinsic = abs(max((sign(q).*(K - F)), 0.0));
    % Map in-the-money to out-of-the-money
    if q * (F - K) > 0
        out = intrinsic + black(F, K, sigma, T, -q);
        return
    end
    out = max(intrinsic, (sqrt(F) * sqrt(K)) * normalised_black(log(F / K), sigma * sqrt(T), q));
end

%% rational cubic related routines

% # Based on
% #
% # â?oShape preserving piecewise rational interpolationâ??, R. Delbourgo, J.A. Gregory - SIAM journal on scientific and
% # statistical computing, 1985 - SIAM. http://dspace.brunel.ac.uk/bitstream/2438/2200/1/TR_10_83.pdf  [caveat emptor:
% # there are some typographical errors in that draft version]
% #
% from __future__ import division
% from math import fabs, sqrt
% 
% from py_lets_be_rational.numba_helper import maybe_jit
% from py_lets_be_rational.constants import DBL_EPSILON
% from py_lets_be_rational.constants import DBL_MIN
% from py_lets_be_rational.constants import DBL_MAX


function out = is_zero(x)
DBL_EPSILON = eps;
DBL_MIN = realmin;
minimum_rational_cubic_control_parameter_value = -(1 - sqrt(DBL_EPSILON));
maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON);

    out = abs(x) < DBL_MIN;
end

function out = rational_cubic_parameter_derivative_left(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l)
    DBL_EPSILON = eps;
    minimum_rational_cubic_control_parameter_value = -(1 - sqrt(DBL_EPSILON));
    h = (x_r - x_l);
    numerator = 0.5 * h * second_derivative_l + (d_r - d_l);
    if is_zero(numerator)
      out = 0;
      return
    end
    denominator = (y_r - y_l) / h - d_l;
    if is_zero(denominator)
      if numerator > 0
      out = maximum_rational_cubic_control_parameter_value;
      else
      out = minimum_rational_cubic_control_parameter_value;
      end
      return
    end
    out = numerator / denominator;
end

  function out = minimum_rational_cubic_control_parameter(d_l, d_r, s, preferShapePreservationOverSmoothness)
    DBL_EPSILON = eps;
    minimum_rational_cubic_control_parameter_value = -(1 - sqrt(DBL_EPSILON));
    monotonic = and( d_l * s >= 0, d_r * s >= 0);
    convex = and(d_l <= s , s <= d_r);
    concave = and(d_l >= s, s >= d_r);
    if and(~monotonic,and(~convex,~concave)) % If 3==r_non_shape_preserving_target, this means revert to standard cubic.
      out = minimum_rational_cubic_control_parameter_value;  
      return 
    end
    d_r_m_d_l = d_r - d_l;
    d_r_m_s = d_r - s;
    s_m_d_l = s - d_l;
    DBL_MAX = realmax;
    r1 = -DBL_MAX;
    r2 = r1;
    maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON);
    % If monotonicity on this interval is possible, set r1 to satisfy the monotonicity condition (3.8).
    if monotonic
        if ~is_zero(s)  % (3.8), avoiding division by zero.
            r1 = (d_r + d_l) / s;  % (3.8)
        elseif  preferShapePreservationOverSmoothness  % If division by zero would occur, and shape preservation is preferred, set value to enforce linear interpolation.
            r1 = maximum_rational_cubic_control_parameter_value;  % This value enforces linear interpolation.
        end
    if or(convex,concave)
        if ~or(is_zero(s_m_d_l),is_zero(d_r_m_s))  % (3.18), avoiding division by zero.
            r2 = max(abs(d_r_m_d_l / d_r_m_s), abs(d_r_m_d_l / s_m_d_l));
        elseif preferShapePreservationOverSmoothness
            r2 = maximum_rational_cubic_control_parameter_value;  % This value enforces linear interpolation.
        end
    elseif and(monotonic,preferShapePreservationOverSmoothness)
        r2 = maximum_rational_cubic_control_parameter_value;  % This enforces linear interpolation along segments that are inconsistent with the slopes on the boundaries, e.g., a perfectly horizontal segment that has negative slopes on either edge.
    end
    out =  max(minimum_rational_cubic_control_parameter_value, max(r1, r2));
    end
  end

function out = rational_cubic_parameter_derivative_right(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r)
    DBL_EPSILON = eps;
    maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON);
    h = (x_r - x_l);
    numerator = 0.5 * h * second_derivative_r + (d_r - d_l);
    if is_zero(numerator)
        out = 0;
        return
    end
    denominator = d_r - (y_r - y_l) / h;
    if is_zero(denominator)
      if numerator > 0
         out = maximum_rational_cubic_control_parameter_value; 
      else
         out = minimum_rational_cubic_control_parameter_value;
         return
      end
    end
    out =  numerator / denominator;
end

function out =  convex_rational_cubic_control_parameter_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r,preferShapePreservationOverSmoothness)
  r = rational_cubic_parameter_derivative_right(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r);
  r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r - y_l) / (x_r - x_l), preferShapePreservationOverSmoothness);
  out =  max(r, r_min);
end

function out = rational_cubic_interpolation(x, x_l, x_r, y_l, y_r, d_l, d_r, r)
    DBL_EPSILON = eps;
    maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON);
    h = (x_r - x_l);
    if abs(h) <= 0
        out =  0.5 * (y_l + y_r);
        return
    end
    % r should be greater than -1. We do not use  assert(r > -1)  here in order to allow values such as NaN to be propagated as they should.
    t = (x - x_l) / h;
    if ~(r >= maximum_rational_cubic_control_parameter_value)
        t = (x - x_l) / h;
        omt = 1 - t;
        t2 = t * t;
        omt2 = omt * omt;
        % Formula (2.4) divided by formula (2.5)
        out = (y_r * t2 * t + (r * y_r - h * d_r) * t2 * omt + (r * y_l + h * d_l) * t * omt2 + y_l * omt2 * omt) / (1 + (r - 3) * t * omt);
        return
    end
    % Linear interpolation without over-or underflow.
    out = y_r * t + y_l * (1 - t);
end

function out =  convex_rational_cubic_control_parameter_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l, preferShapePreservationOverSmoothness)
  r = rational_cubic_parameter_derivative_left(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l);
  r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r - y_l) / (x_r - x_l), preferShapePreservationOverSmoothness);
  out = max(r, r_min);
end


