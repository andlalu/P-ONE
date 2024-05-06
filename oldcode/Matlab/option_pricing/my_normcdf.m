function p = my_normcdf(x)
% DESCRIPTION O BE ADDED AT A LATER STAGE
% In a nutshell, faster version of normcdf available in matlab. All copyrights are as in the original MATLAB code

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p = 0.5 * erfc(-x ./ sqrt(2));

end
