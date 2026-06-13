%%%%% 11.08.2014 v.1.0 - Beta for: SVHJ; Alpha for: any other model;

function [states,x] = impl_state_func_fut(deriv_data, TT_marks, type, theta, varargin)
%IMPL_STATE_FUNC Implies the level of the unobservable states from a set of market implied volatilities and futures prices
%   Function finds the level of unobservable states which minimize Black-Scholes implied volatility pricing errors 
%       between the theoretical option price given the parameters and the observed market prices 
%   INPUTS:
%   deriv_data = matrix of derivatives prices, risk free and dividend rates with the following structure:
%  'deriv_data' structure: ||1 sample period|2 effective maturity (in years)|3 strike_price (absolute)|4 implied vol |5 futures price |6 interest rate||

%% This part sets-up the state backing out according to the model type:
switch type 
    case 'SV'
        NS=1;
        states_start=[theta(3)];
        states_min=[0.0001];
        states_max=[1];
    case 'SVJ'
        NS=1;
        states_start=[theta(6)];
        states_min=[0.0001];
        states_max=[1];
    case 'SVHJ'
        NS=2; %number of states 
        states_start=[theta(6),theta(10)];
        states_min=[0.0001,theta(10)];
        states_max=[.3,250];
    case 'SVVJ' 
        NS=1;
        states_start=[theta(6)];
        states_min=[0.0001];
        states_max=[1];
end

%% Compute the characteristic function coefficients for all the maturity points of this sample
TT=unique(deriv_data(:,2));
M=max(deriv_data(:,1));
if nargin==5
    cf_coeff=varargin{1};
else
cf_coeff = price_fft_fut(0,0,TT,0,type,theta,0,1);
end
%% Implying states for each day and each maturity in the sample
options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6); %Optimize silently, w/ no display in the command window %Optional:'TolFun',1e-12,'TolX',1e-8
states=zeros([M,NS]);
x=zeros([M,1]);
%Imply states on 1st date starting from long run averages
aux=deriv_data(deriv_data(:,1)==1,:);
func = @(state) vector_errors(state, aux, TT_marks(1,2:end),type,theta,cf_coeff);  
[states(1,:),~,~,x(1,:)]=lsqnonlin(func,states_start,states_min,states_max,options);

for i=2:M
aux=deriv_data(deriv_data(:,1)==i,:);
func = @(state) vector_errors(state, aux, TT_marks(1,2:end),type,theta,cf_coeff);  
[states(i,:),~,~,x(i,:)]=lsqnonlin(func,states(i-1,:),states_min,states_max,options);
end
states=[diff(log(deriv_data([1;1+cumsum(TT_marks(1:M-1,end))],5))),states(2:end,:)];
end

%% Auxilliary function that computes the pricing error differences
function out = vector_errors (state,aux,TT_mark,type,theta,cf_coeff)
out=zeros([size(aux(:,1),1),1]);
%     aux2=price_fft_fut([log(aux(TT_mark(t),5));state';], aux(TT_mark(t-1)+1:TT_mark(t),3).*aux(TT_mark(t),5), aux(TT_mark(t),2),aux(TT_mark(t),6),type,theta,1,0,cf_coeff(:,:,t-1));

for t=2:size(TT_mark,2)
    aux2=price_fft_fut([log(aux(TT_mark(t),5));state';], aux(TT_mark(t-1)+1:TT_mark(t),3).*aux(TT_mark(t),5), aux(TT_mark(t),2),aux(TT_mark(t),6),type,theta,1,0,cf_coeff(:,:,t-1));
    out(TT_mark(t-1)+1:TT_mark(t))=aux2(:,2);
end
aux2=aux(:,4);
out=out-aux2;
end

