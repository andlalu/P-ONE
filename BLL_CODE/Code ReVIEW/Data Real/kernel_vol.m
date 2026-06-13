function [ ivols ] = kernel_vol( gridpoints, deriv_data, varargin )
%kernel_vol computes a (Gaussian) kernel based estimator for the the implied volatilty at a certain maturity and moneyness level provided the estimator 
%   INPUT
%   gridpoints - array containing the maturities and moneyness levels for which the kernel estimated implied volatilites should be provided 
%              - structure: |1 maturity (in days)|2 strike level 
%   deriv_data - array containing maturities, moneyness levels and implied volatilities accross all options of that specific day
%   varargin   - possibly the current stock price to compute the moneyness level
%  
%   OUTPUT
%   ivols - column vector containing the kernel-estimated implied volatilities at the maturity/moneyness levels contained in gridpoints

%Kernel weights set-up
h_T= 0.01; %kernel weight for maturity (T)
h_m= 0.001; %kernel weight for moneyness 

N=length(gridpoints); %the number of output points
D=length(deriv_data); %the number of derivatives used as input

gridpoints=reshape(repmat(reshape(gridpoints,[1,2*N]),[D,1]),[N*D,2]);
deriv_data=repmat(deriv_data,[N,1]);
kernel1=exp(-((gridpoints(:,1)-deriv_data(:,1)).^2)./(2*h_T));
kernel2=exp(-((gridpoints(:,2)-deriv_data(:,2)).^2)./(2*h_m));
kernel_den=reshape(kernel1.*kernel2,[D,N]); 
kernel_num=reshape(deriv_data(:,3).*kernel1.*kernel2,[D,N]);
ivols=(sum(kernel_num,1)./sum(kernel_den,1))';

end

