%%%% 12.08.2014 V2.0  Beta for: any model;
function [ answ ] = estimate_hess2( deriv_data, TT_marks, dt, type, theta,lb,ub,theta1step,X )    
    

  %Criterion characteristics
  factor=1e3;
  alpha=0.01;
  
  %Minimization routine
  func = @(theta) aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha,theta1step);
%   answ=fminsearch(func,boundtofree(theta,lb,ub),optimset('Display','iter','TolFun',1e0,'MaxFunEvals',1000));
%   answ = freetobound(answ,lb,ub)
  [answ,answ2]=func(boundtofree(type,theta,lb,ub));
end

function [out, X] = aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha,theta1step)
theta=freetobound(type,theta,lb,ub);
disp(theta);
disp(theta1step);
X=impl_state_func_fut(deriv_data, TT_marks, type,theta); 
H0=funcH0(theta1step, X, dt, type, factor);
out=critstep1(theta,H0, X, dt, type,1,factor,alpha)
end