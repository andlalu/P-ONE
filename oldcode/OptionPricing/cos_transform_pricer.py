import numpy as np

""" Option pricing routine"""


def cos_transform_pricer():
	a = -5.0  # integration lower bound
	b = 5.0  # integration upper bound
	c = .0  # other integration bound
	d = b  # integration bound
	n_cos = 1024  # Integration bound	for COS pricing part
	k = np.linspace(0, n_cos-1, n_cos)  # integration grid values
	# Routine set - up:
	arg_dma = k * np.pi * (d - a) / (b - a)
	arg_cma = k * np.pi * (c - a) / (b - a)
	fact_sin = k * np.pi / (b - a)
	chi = np.cos(arg_dma) * np.exp(d) - np.cos(arg_cma) * np.exp(c)
	chi = chi + np.multiply(fact_sin, np.sin(arg_dma)) * np.exp(d)
	chi = chi - np.multiply(fact_sin, np.sin(arg_cma)) * np.exp(c)
	test = np.power((1 + np.square(fact_sin)), (-1))
	chi = np.multiply(chi, test)
	phi = np.multiply((np.sin(arg_dma) - np.sin(arg_cma)),
							 (np.append(np.multiply((b - a), np.power((k[1:] * np.pi), (-1))), 0) +
								np.append(np.zeros((n_cos - 1, 1)), d - c)))
	v_call = 2.0 / (b - a) * (chi - phi)
	return v_call # assume okay for now

# u = ((cumsum(ones(N_cos, 1)) - 1) * pi / (b - a)) * 1i; # integration grid	coordinates

# 	if (cf_coeff_save == 1)
#     ## Function returns ccf coefficients
#     out = nan([N_cos, len, length(T)]);
#     parfor
#     i = 1:N_cos
#     out(i,:,:) = char_func(k(i) * pi / (b - a));
#
# 	else if(cf_coeff_save == 0)
# 	## Find option prices
# 	cf_coeff = varargin{1};
# 	out = nan([length(K), length(T) + 1]);
# 	out(:, 1) = K;
# for tt=1:length(T)
# for kk=1:length(K)
# cf = exp(cf_coeff(:, 1, tt) + cf_coeff(:, 2: end, tt)*X0 - 1
# i * k * pi / (b - a) * (log(K(kk))));
# out(kk, tt + 1) = exp(-r(tt) * T(tt)) * K(kk) * real(sum([.5;
# ones(N_cos - 1, 1)].*cf. * exp(-1
# i * (k * a * pi / (b - a))).*V_call));
# end
# end
# % Return
# prices in IV / "Dollar"
# if out_BS_vol == 1
#     for tt=1:length(T)
#     for kk=1:length(K)
#     out(kk, tt + 1) = bs_imp_vol_fut(abs(out(kk, tt + 1)), T(tt), exp(X0(1)), K(kk), r(tt));
# end
# end
#
# end
# end
