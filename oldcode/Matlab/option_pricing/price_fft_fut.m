    %%%%%% v 2.0 Prices options via FFT for a set of AJD types:

function out = price_fft_fut(X0, K, T, r, type, theta, out_BS_vol, cf_coeff_save, varargin)
%CHAR_FUNC Price call options via conditional characteristic functions using Fast Fourier Transform algorithms (Carr&Madan JCF 1999)
%   INPUTS:
%   X0= vector of initial states: e.g. for Hawkes processes (\logS_0, v_0, \lambda_0)
%   K = vector of strike prices for which option prices should be determined (can be of any arbitrary length)
%   T = scalar indicating the end/final time at which the transform is evaluated
%   type = string, denoting model type, e.g. 'SVHJ'
%   theta = model parameters
%   out_BS_vol = '1' - provide output prices expressed in B-S implied volatilities instead of dollar-value prices
%   cf_coeff_save = '1' - save the CF coefficients and provide them as an output variable instead of the price ranges.

%   OPTIONAL INPUTS:
%   cf_coeff_inp = Matrix of cf_coefficients used to price the option instead of computing them within the function 

%NB:  THIS FUNCTION IS OPTIMIZED FOR PARALLEL & SINGLE THREAD EXECUTION

%% This part selects the characteristic function according to the model type: 
switch type 
    %%%% UNIVARIATE MODEL SPECIFICATIONS.
    case 'SVJ' 
        len=3;
        % to be added. 
         %Allocate parameters
        muj=theta(1);
        mujq=theta(2);
        sigmaj=theta(3);
        eta=theta(4);
        k=theta(5);
        vbar=theta(6);
        sigma=theta(7);
        rho=theta(8);
        lambda=theta(9);
        %Define CCF system matrixes
        K0 = [ -(exp(mujq+0.5*sigmaj^2)-1)*lambda;      k*vbar;];
        K1 = [ 0, -0.5;
            0,  -k;];
        H00 = [ 0,      0;
            0,      0;];
        H1(:,:,1) = [ 0,0; 0,0;];
        H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
        L0=lambda;
        L1=[0;0;];
        JT=@(beta) (exp(beta(1)*mujq+0.5*beta(1)^2*sigmaj^2)-1);
        %Allocate function handle
        char_func = @(arg) ccf([arg;0],T,K0,K1,H00,H1,L0,L1,JT);
        
    case 'SVHJ'
        len=4;
        %Allocate parameters    
        muj=theta(1);
        mujq=theta(2);
        sigmaj=theta(3);
        eta=theta(4);
        k=theta(5);
        vbar=theta(6);
        sigma=theta(7);
        rho=theta(8);
        kl=theta(9);
        lambda=theta(10);
        delta=theta(11);
        %Define CCF system matrixes
        c=1;
        K0  = [ 0;    k*vbar;   kl*lambda];
        K1  = [ 0,   -0.5,  -c*(exp(mujq+0.5*sigmaj^2)-1);
                0,     -k,                            0;
                0,      0,                          -kl; ];
        H00 = [ 0,   0,   0;
                0,   0,   0;
                0,   0,   0; ];
        H1(:,:,1) = [ 0,0,0; 0,0,0; 0,0,0;];
        H1(:,:,2) = [ 1,sigma*rho, 0; sigma*rho,sigma^2,0; 0,0,0;];
        H1(:,:,3) = [ 0,0,0; 0,0,0; 0,0,0;];
        L0 = 0;
        L1 = [0;0;c];
        JT = @(beta) (exp(beta(1)*mujq+0.5*beta(1).^2*sigmaj^2+delta*beta(3))-1);
        %Allocate function handle
        char_func = @(arg) ccf([arg;0;0],T,K0,K1,H00,H1,L0,L1,JT);
        
        %%%%% MULTIVARIATE SPECIFICATIONS
        case '2SVHJ'
        len=5;
        %Allocate parameters
        [muj,mujq,sigmaj,eta,k,vbar,sigma,rho,kl1,lambda1,delta1,mut_delta1,kl2,lambda2,delta2,mut_delta2] = ...
            deal(theta(:,1),theta(:,2),theta(:,3),theta(:,4),theta(:,5),theta(:,6),...
            theta(:,7),theta(:,8),theta(:,9),theta(:,10),theta(:,11),theta(:,12),theta(:,13),theta(:,14),theta(:,15),theta(:,16));
        %Define CCF system matrixes
        c=1;
        K0  = [ 0; k.*vbar;   kl1.*lambda1; kl2.*lambda2;];
        K1  = [ 0,-0.5, -c*(exp(mujq(1)+0.5*sigmaj(1).^2)-1),  0;
                0,-k(1), 0,                                    0;  
                0, 0,   -kl1,   0; 
                0, 0,    0,      -kl2; ];
            
        H00 = zeros(4,4);
        H1(:,:,1) = [ 0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0;];
        H1(:,:,2) = [ 1, sigma*rho, 0, 0; sigma*rho, sigma^2, 0, 0; 0,0,0 0; 0,0,0 0;];
        H1(:,:,3) = [ 0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0; ];
        H1(:,:,4) = [ 0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0; ];
        L0 = [0 0];
        L1 = [0 0; 0 0; c 0; 0 c;];
        JT = @(beta) ([exp(beta(1)*mujq(1)+0.5*beta(1).^2*sigmaj(1)^2+delta1*beta(3)+mut_delta1*beta(4))-1; exp(delta2*beta(4)+mut_delta2*beta(3))-1]);
        %Allocate function handle
        char_func = @(arg) ccf([arg;0;0;0;],T,K0,K1,H00,H1,L0,L1,JT);
        
        case '2SVhatHJ'
        len=4;
        %Allocate parameters
        [muj,mujq,sigmaj,eta,kl1,lambda1,delta1,mut_delta1,kl2,lambda2,delta2,mut_delta2,vhat] = ...
            deal(theta(1),theta(2),theta(3),theta(4),theta(5),theta(6),theta(7),theta(8),theta(9),theta(10),theta(11),theta(12),theta(13));
        %Define CCF system matrixes
        c=1;
        K0  = [ -vhat/2; kl1.*lambda1; kl2.*lambda2;];
        K1  = [ 0,-c*(exp(mujq(1)+0.5*sigmaj(1).^2)-1),  0;  
                0, -kl1,   0; 
                0, 0,      -kl2; ];
            
        H00 = [vhat 0 0; 0 0 0; 0 0 0;];
        H1(:,:,1) = [ 0,0,0; 0,0,0; 0,0,0;];
        H1(:,:,2) = [ 0,0,0; 0,0,0; 0,0,0;];
        H1(:,:,3) = [ 0,0,0; 0,0,0; 0,0,0;];
        L0 = [0 0];
        L1 = [0 0; c 0; 0 c;];
        JT = @(beta) ([exp(beta(1)*mujq(1)+0.5*beta(1).^2*sigmaj(1)^2+delta1*beta(2)+mut_delta1*beta(3))-1; exp(delta2*beta(3)+mut_delta2*beta(2))-1]);
        %Allocate function handle
        char_func = @(arg) ccf([arg;0;0;],T,K0,K1,H00,H1,L0,L1,JT);
end

%% This part of the function sets the dampening factor and integration grid 'roughness'
alpha = 1.2;
eta = 0.1; 
N = 2^13; %Must be a power of 2, e.g. 4096=2^12; !!! IN BETA CANNOT CHANGE THIS ONLY IN THIS FUNCTION!!!!!!!!!

%% Computing Option Prices using FFT
%Construct the integrand function & use Simpson's rule to determine integration points
v = 0:eta:(N-1)*eta;
lambda = 2*pi/(N*eta);
k = -(N*lambda/2):lambda:(-N*lambda/2+(N-1)*lambda);

%Construct integrand for the FFT pricing of in-the-money options
if nargin<9 %check if CF function coefficients aren't provided by the user
    if cf_coeff_save  %check if instead of prices, the function is requested to output the cf coefficient matrix
        out=zeros([N,len,length(T)]);
        parfor i=1:N
            out(i,:,:) = char_func(v(i)-(1+alpha)*1i);
        end   
    else
        error('This set-up of the function was not programmed'); %This would be inefficient, so this set-up is not programmed.
    end
else %if the cf function coefficients were provided then use them to build the integrand and evaluate the derivative prices 
    cf_coeff=varargin{1};
        weights= eta/3.*(3.*ones([N,1])+((-1)*ones([N,1])).^(cumsum(ones([N,1])))-[1;zeros([N-1,1])]); 
        integrand = repmat(exp(-r.*T'),[N,1]).*permute(repmat(exp(1i*N*lambda/2.*v')./((alpha+1i*v').*(1+alpha+1i*v')),[1,1,length(T)]).*exp(cf_coeff(:,1,:)+sum(cf_coeff(:,2:len,:).*repmat(X0',[N,1,length(T)]),2)).*repmat(weights,[1,1,length(T)]),[1,3,2]);
    aux(:,1)=exp(k)';
    aux(:,2:2+length(T)-1)=real(repmat(exp(-alpha*k')/pi,[1,length(T)]).*fft(integrand));
    out(:,1)=K;
    out(:,2:2+length(T)-1)=(interp1q(aux(:,1),aux(:,2:end),K));   
    if out_BS_vol %check if B-S implied volatilities should be outputted instead of dollar value prices
            aux_S=exp(X0(1,1));
            for t=2:length(T)+1
            for i=1:size(out,1)
               out(i,t) = bs_imp_vol_fut(out(i,t), T(t-1), aux_S, K(i),r); 
            end
            end
    end
end

end


