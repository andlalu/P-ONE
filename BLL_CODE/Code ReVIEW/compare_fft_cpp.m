%C++ code comparison with matlab. 
  muj0=-0.05;
  mujq0=-0.1;
  sigmaj0=0.2;
  eta0=2;
  k0=5;
  vbar0=0.01;
  sigma0=0.25;
  rho0=-0.5;
  kl0=15; %Starting value is the difference from delta
  lambda0=0.5;
  delta0=10;
  theta = [muj0,mujq0,sigmaj0,eta0,k0,vbar0,sigma0,rho0,kl0,lambda0,delta0];%theta SVHj
  
  X0 = [log(100);0.04;10];
  K = [.95; 1; 1.05]*exp(X0(1));
  T = [0.1; 0.5; 0.95];
  r = zeros(3,1);
  type= 'SVHJ';
  coeff_save = price_fft_fut(X0, K, T, r, type, theta, 0, 1);
  prices = price_fft_fut(X0, K, T, r, type, theta, 0, 0, coeff_save);