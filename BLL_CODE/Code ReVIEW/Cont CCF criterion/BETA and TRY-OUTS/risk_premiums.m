%% Computing fancy risk premium based on the model: 
clearvars vrp erp_tot erp_j erp_v
for i=1:length(states)
T=22/252;
erp_v(i,1) =eta0*(vbar0+(((1-exp(-k0*T))*(states(i,2)-vbar0))/(k0*T)));

erp_j(i,1) =(exp(muj0+0.5*sigmaj0^2)-exp(mujq0+0.5*sigmaj0^2))*( (kl0*lambda0)/(kl0-delta0) + (1-exp(-(kl0-delta0)*T))*(((kl0-delta0)*states(i,3)-kl0*lambda0)/(((kl0-delta0).^2)*T)));

erp_tot=erp_j+erp_v;

vrp(i,1) =(muj0^2-mujq0^2)*( (kl0*lambda0)/(kl0-delta0) + (1-exp(-(kl0-delta0)*T))*(((kl0-delta0)*states(i,3)-kl0*lambda0)/(((kl0-delta0).^2)*T)));
end
ma_mat = repmat(1:length(states)-3,[4,1])+repmat((0:3)',[1,length(states)-3]);
ma_erp = sum(erp(ma_mat),1)./4;
ma_vrp = sum(vrp(ma_mat),1)./4;
