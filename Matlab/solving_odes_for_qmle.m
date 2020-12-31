clc
clear
syms etaSTAR vzero kappa vbar sigmav rho u(t)
eqn = diff(u(t),t) == (etaSTAR) * (vzero^2*exp(-2*kappa*t) + (2*kappa*vbar + sigmav^2)/(kappa)*(exp(-kappa*t)*(vzero-vbar)*(1-exp(-kappa*t)) + vbar/2*(1-exp(-2*kappa*t))))...
                      + kappa*vbar*((etaSTAR)*(vbar*t + (vzero-vbar)*(1/kappa)*(1-exp(-kappa*t))))...
                      + sigmav*rho*(exp(-kappa*t)*(vzero-vbar)+vbar)...
                      - kappa*u(t);
cond = [u(0)==0];
EYV(t) = dsolve(eqn,cond);
% (exp(-2*kappa*t)*(exp(2*kappa*t)*((etaSTAR*sigmav^2*vbar)/kappa + 2*rho*sigmav*vbar - 2*etaSTAR*vbar^2 + 2*etaSTAR*vzero*vbar) - 2*etaSTAR*vbar^2 - 2*etaSTAR*vzero^2 - t*exp(kappa*t)*(2*etaSTAR*sigmav^2*vbar - 2*etaSTAR*vzero*sigmav^2 + 2*kappa*rho*sigmav*vbar - 2*kappa*rho*vzero*sigmav + 2*etaSTAR*kappa*vbar^2 - 2*etaSTAR*kappa*vzero*vbar) + 4*etaSTAR*vbar*vzero - (etaSTAR*sigmav^2*vbar)/kappa + (2*etaSTAR*sigmav^2*vzero)/kappa + 2*etaSTAR*kappa*t*vbar^2*exp(2*kappa*t)))/(2*kappa) - exp(-kappa*t)*(((2*etaSTAR*sigmav^2*vzero)/kappa + 2*rho*sigmav*vbar - 4*etaSTAR*vbar^2 + 6*etaSTAR*vbar*vzero - 2*etaSTAR*vzero^2)/(2*kappa) - vzero^2)
EY = etaSTAR*(vbar*t + (vzero-vbar)/kappa*(1-exp(-kappa*t)));
EV = exp(-kappa*t)*(vzero-vbar)+vbar;
EV2 = vzero^2 * exp(-2*kappa*t) + (2*kappa*vbar + sigmav^2)/(kappa) * (exp(-kappa*t)*(vzero-vbar)*(1-exp(-kappa*t))+vbar/2*(1-exp(-2*kappa*t)));
VAR_V = simplify(EV2 - power(EV,2));
COV_YV = simplify(EYV-EY*EV);
% VAR_Y = 
eqn = diff(u(t),t) == 2*etaSTAR*EYV + EV;
cond = [u(0) == 0]; %> is this right?
EY2(t) = simplify(dsolve(eqn,cond));
VAR_Y = simplify(EY2 - power(EY,2)); 
