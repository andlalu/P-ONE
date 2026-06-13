%%%% 12.08.2014 V2.0  Beta for: any model;
function [ answ ] = estimate_hess( deriv_data, TT_marks, dt, type, theta,lb,ub,theta1step,X,H0)    
    

  %Criterion characteristics
  factor=1e3;
  alpha=0.01;
  
  %Minimization routine
  func = @(theta) aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha,X,H0);
%   answ=fminsearch(func,boundtofree(theta,lb,ub),optimset('Display','iter','TolFun',1e0,'MaxFunEvals',1000));
%   answ = freetobound(answ,lb,ub)
  [answ,answ2]=func(boundtofree(type,theta,lb,ub));
end

function [out, X] = aux(deriv_data, TT_marks, dt, type, theta,lb,ub,factor,alpha,X,H0)
theta=freetobound(type,theta,lb,ub);
disp(theta);
out=critstep1(theta,H0, X, dt, type,1,factor,alpha)
end