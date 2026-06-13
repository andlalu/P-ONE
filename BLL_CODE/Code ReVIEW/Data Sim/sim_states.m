%%%% 01.06.2015 V2.0 Final for: SV,SVJ,SVHJ; Alpha for: COX, SVVJ;

function [ states,varargout ] = sim_states(X0,N,dt,type,theta,varargin)
%SIMSERIES 
%         Optional inputs:
%           M = number of random draws in between two data points, default number is M=100
%           seed = the seed of the random number generator 

    M=500; % number of random draws in between two data-points
    burnin=0;
    P=M*(N+burnin); % total number of points drawn along the simulated path of the state series
    deltat=dt/M;
if max(size(varargin))>0
    seed=varargin{1};
else
    seed=1+floor(1e6*rand(1,1));
end
f = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(f);
fprintf('Random generator seed: %u\n',seed);

switch type
    case 'SVHJ'
    %Allocate parameters
    muj=theta(1);
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigmav=theta(7);
    rho=theta(8);
    kl=theta(9);
    lambdabar=theta(10);
    delta=theta(11);
    %Initialize the vectors which will store the states
    s=[X0(1); zeros([P,1])]; 
    v = [X0(2); zeros([P,1])]; 
    lambda=[X0(3); zeros([P,1])];

    %Draw the Brownian components and the jumps "instantaneous CDF"
    W = sqrt(deltat).*randn([P,2]);
    U = rand([P,1]); 

    %Calculate the path
    for i=2:P+1 
    v(i) = max(v(i-1)+k*(vbar-v(i-1))*deltat + sigmav*sqrt(v(i-1))*(rho*W(i-1,1)+sqrt(1-rho^2)*W(i-1,2)),0);%Discretized SVol WITH ABSORPTION!
    lambda(i) = lambda(i-1)+ kl*(lambdabar-lambda(i-1))*deltat + delta*(U(i-1)<lambda(i-1)*deltat);
    s(i) = s(i-1) +((-exp(mujq+0.5*sigmaj^2)+1)*lambda(i-1)+(eta-0.5)*v(i-1))*deltat + sqrt(v(i-1))*W(i-1,1) + (U(i-1)<lambda(i-1)*deltat)*(muj+sigmaj.*randn([1,1])); %Discretized log stock path
    end
    s=s(burnin*M+1:M:P);
    v=v(burnin*M+1:M:P);
    lambda=lambda(burnin*M+1:M:P);

    states=[s,v,lambda];
    varargout{1}=seed;
    
    case 'SVJ'
    %Allocate parameters: 
    muj=theta(1);
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigmav=theta(7);
    rho=theta(8);
    lambda=theta(9);
    %Initialize the vectors which will store the states
    s = [X0(1); zeros([P,1])]; 
    v = [X0(2); zeros([P,1])]; 

    %Draw the Brownian components and the jumps "instantaneous CDF"
    W = sqrt(deltat).*randn([P,2]);
    U = rand([P,1]); 

    %Calculate the path
    for i=2:P+1 
    v(i) = max(v(i-1)+k*(vbar-v(i-1))*deltat + sigmav*sqrt(v(i-1))*(rho*W(i-1,1)+sqrt(1-rho^2)*W(i-1,2)),0);%Discretized SVol WITH ABSORPTION!
    s(i) = s(i-1) + ((-exp(mujq+0.5*sigmaj^2)+1)*lambda+(eta-0.5)*v(i-1))*deltat + sqrt(v(i-1))*W(i-1,1) + (U(i-1)<lambda*deltat)*(muj+sigmaj.*randn([1,1])); %Discretized log stock path
    end
    states=[s,v];
    varargout{1}=seed;
    
    case 'SV'
    %Allocate parameters
    eta=theta(1);
    k=theta(2);
    vbar=theta(3);
    sigmav=theta(4);
    rho=theta(5);
    %Initialize the vectors which will store the states
    s = [X0(1); zeros([P,1])]; 
    v = [X0(2); zeros([P,1])]; 

    %Draw the Brownian components and the jumps "instantaneous CDF"
    W = sqrt(deltat).*randn([P,2]); 

    %Calculate the path
    for i=2:P+1 
    v(i) = max(v(i-1)+k*(vbar-v(i-1))*deltat + sigmav*sqrt(v(i-1))*(rho*W(i-1,1)+sqrt(1-rho^2)*W(i-1,2)),0);%Discretized SVol WITH ABSORPTION!
    s(i) = s(i-1) - 0.5*v(i-1)*deltat + eta*v(i-1)*deltat + sqrt(v(i-1))*W(i-1,1); %Discretized stock path
    end
    varargout{1}=seed;
    s=s(burnin*M+1:M:P);
    v=v(burnin*M+1:M:P);
    states=[s,v];
    
    case 'PAN'
        warning('Not tested yet');
    %Allocate parameters
    muj=theta(1);
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigmav=theta(7);
    rho=theta(8);
    lambda0=theta(9);
    lambda1=theta(10);
    %Initialize the vectors which will store the states
    s=[X0(1); zeros([P,1])]; 
    v = [X0(2); zeros([P,1])]; 
    lambda=[lambda0+lambda1*X0(2);zeros([P,1])];

    %Draw the Brownian components and the jumps "instantaneous CDF"
    W = sqrt(deltat).*randn([P,2]);
    U = rand([P,1]); 

    %Calculate the path
    for i=2:P+1 
    v(i) = max(v(i-1)+k*(vbar-v(i-1))*deltat + sigmav*sqrt(v(i-1))*(rho*W(i-1,1)+sqrt(1-rho^2)*W(i-1,2)),0);%Discretized SVol WITH ABSORPTION!
    lambda(i) = lambda0+lambda1*v(i);
    s(i) = s(i-1) +((-exp(mujq+0.5*sigmaj^2)+1)*lambda(i-1)*deltat+(eta-0.5)*v(i-1))*deltat + sqrt(v(i-1))*W(i-1,1) + (U(i-1)<lambda(i-1)*deltat)*(muj+sigmaj.*randn([1,1])); %Discretized log stock path
    end
    s=s(burnin*M+1:M:P);
    v=v(burnin*M+1:M:P);
    lambda=lambda(burnin*M+1:M:P);
%     varargout{1}=lambda;
    states=[s,v,lambda];
    varargout{1}=seed;    
    
    case 'COX'
    error('Not tested yet');
    %Allocate parameters
    muj=theta(1);
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigmav=theta(7);
    rho=theta(8);
    kl=theta(9);
    lambdabar=theta(10);
    delta=theta(11);
    %Initialize the vectors which will store the states
    s=[X0(1); zeros([P,1])]; 
    v = [X0(2); zeros([P,1])]; 
    lambda=[X0(3); zeros([P,1])];

    %Draw the Brownian components and the jumps "instantaneous CDF"
    W = sqrt(deltat).*randn([P,2]);
    U = rand([P,1]); 

    %Calculate the path
    for i=2:P+1 
    v(i) = max(v(i-1)+k*(vbar-v(i-1))*deltat + sigmav*sqrt(v(i-1))*(rho*W(i-1,1)+sqrt(1-rho^2)*W(i-1,2)),0);%Discretized SVol WITH ABSORPTION!
    lambda(i) = kl*(lambdabar-lambda(i-1))*deltat + delta*(U(i-1)<lambda(i-1)*deltat);
    s(i) = s(i-1) +((-exp(mujq+0.5*sigmaj^2)+1)*lambda(i-1)+(eta-0.5)*v(i-1))*deltat + sqrt(v(i-1))*W(i-1,1) + (U(i-1)<lambda(i-1)*deltat)*(muj+sigmaj.*randn([1,1])); %Discretized log stock path
    end
    s=s(burnin*M+1:M:P);
    v=v(burnin*M+1:M:P);
    lambda=lambda(burnin*M+1:M:P);
%     varargout{1}=lambda;
    states=[s,v,lambda];
    varargout{1}=seed;  
end

end


