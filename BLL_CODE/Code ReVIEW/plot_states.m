figure
set(0,'defaulttextinterpreter','latex');

%Plot first state
subplot(3,1,1);
plot(1:500,100*diff(states(1:end,1)));
legend('$\log\left(y_{t+\Delta}\right)-\log\left(y_t\right)$','Interpreter','latex');
xlabel('Observation (week) no.');
ylabel('\%');
title('Simulated state plots:');
hold on
%Plot second state
subplot(3,1,2);
plot(1:500,100*states(2:end,2));
legend('$v$');
xlabel('Observation (week) no.');
ylabel('\%');
hold on 
%Plot second state
subplot(3,1,3);
plot(1:500,states(2:end,3));
legend('$\lambda$');
xlabel('Observation (week) no.');
hold off