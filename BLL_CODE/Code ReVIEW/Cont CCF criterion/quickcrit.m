%%%% 12.08.2014 V2.0  Beta for: any model;
function [ answ ] = quickcrit( deriv_data, TT_marks, dt, type, theta,lb,ub )    
    

  %Criterion characteristics
  factor=10;
  alpha=0.01;
  
  %Minimization routine first step
  func = @(theta) aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha);
%   answ=fminsearch(func,boundtofree(type,theta,lb,ub),optimset('Display','iter','TolFun',1e0,'MaxFunEvals',1000));
%   answ = freetobound(type,answ,lb,ub)
    answ=func(theta);
end

function [out] = aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha)
% theta=freetobound(type,theta,lb,ub);
disp(theta);
X=impl_state_func_fut(deriv_data, TT_marks, type,theta); 
H0=funcH0(theta, X, dt, type, factor);
out=critstep1(theta,H0, X, dt, type,0,factor,alpha)
end