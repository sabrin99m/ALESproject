clc
clearvars
close all

rng('default')

% Simulation parameters

N = 1000;   % NUmber of data 
L = 2;      % Number of coefficients of the system

% True system to be estimated
theta1_true = ones(N,1);    

for i = 401:N
    theta1_true(i) = sin(0.01*pi*i + pi/2);
end

theta2_true = ones(N,1);

sigma_u = 10;
u1 = sigma_u*randn(N,1);    
u2 = sigma_u*randn(N,1);
u = [u1 u2];
sigma_v = 1;
v = sigma_v*randn(N,1);        

y_noiseless = theta1_true.*u1 + theta2_true.*u2;
y = theta1_true.*u1 + theta2_true.*u2 + v;  % 'Unknown' system output

% Four estimation methods are compared: 
% plain RLS (lambda = 1), 
% RLS with constand Forgetting Factor (lambda = 0.95)
% GVFF RLS
% Variable Forgetting Factor RLS 

% Plain RLS (form 3) 
P = 10^6*eye(L);
theta_rls = zeros(N,L)';

for i = L:N
    b = 1 + u(i,:)*P*u(i,:)';
    P = P - (1/b)*P*u(i,:)'*u(i,:)*P;
    K = P*u(i,:)';
    e = y(i) - u(i,:)*theta_rls(:,i-1);
    theta_rls(:,i) = theta_rls(:,i-1) + K*e;
end

theta_rls = theta_rls';

% RLS with constant forgetting factor (lambda = 0.95)

P = 10^6*eye(L);
theta_crls = zeros(N,L)';
lambda = .95;
for i = L:N
    b = lambda + u(i,:)*P*u(i,:)';
    P = (1/lambda)*(P - (1/b)*P*u(i,:)'*u(i,:)*P);
    K = P*u(i,:)';
    e = y(i) - u(i,:)*theta_crls(:,i-1);
    theta_crls(:,i) = theta_crls(:,i-1) + K*e;
end

theta_crls = theta_crls';


% GVFF
%___________________



% Variable Forgetting Factor RLS 

P = 10^6*eye(L);
theta_vff = zeros(N,L)';
lambda_max = 0.999999;                 % Maximum forgetting factor
gamma = 1.5;                           % lamba(n) = lambda_max when sigma_e(n) <= gamma*sigma_v(n)
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates

K_alpha = 2;
K_beta = 5*K_alpha;
alpha = 1 - 1/(K_alpha*L);
beta = 1 - 1/(K_beta*L);

% Initialization of standard deviations
q = u(2,:) * P *u(2,:)';               
sigma_q = sqrt(var(q));                % prediction error associated with the input u(n)
sigma_e = sqrt(mean(y_noiseless.^2));  % power of the a priori error signal

lambda_n = ones(N, 1) * lambda_max;    % Variable Forgetting Factor initialization

for i = 2:N
    e = y(i) - u(i,:)*theta_vff(:,i-1);                            % A priori error
    K = (P * u(i,:)') / (lambda_n(i-1) + u(i,:) * P * u(i,:)');  % Kalman gain vector
    theta_vff(:,i) = theta_vff(:,i-1) + K*e;
    q_n = u(i,:) * P *u(i,:)';
    P = (1/lambda_n(i))*(P - K*u(i,:)*P);

    % Update forgetting factor
    sigma2_e = alpha*sigma_e^2 + (1-alpha)*e^2;
    sigma_e = sqrt(sigma2_e);
    sigma2_q = alpha*sigma_q^2 + (1-alpha)*q_n^2;
    sigma_q = sqrt(sigma2_q);
    
    sigma2_v = beta*sigma_v^2 + (1-beta)*e^2;
    sigma_v = sqrt(sigma2_v);

    lambda_n(i) = min((sigma_q*sigma_v)/(xi + abs(sigma_e - sigma_v)), lambda_max);
end

theta_vff = theta_vff';


% Plot results
figure
title('Plain RLS vs. VFF-RLS')
t = 1:N;    % temporal line for plot

subplot(3,2,1)
plot(t,theta1_true)
hold on 
plot(t, theta_rls(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'plainRLS')

subplot(3,2,2)
plot(t,theta2_true)
hold on 
plot(t, theta_rls(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'plainRLS')

subplot(3,2,3)
plot(t,theta1_true)
hold on 
plot(t, theta_crls(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLS with mu = 0.95')

subplot(3,2,4)
plot(t,theta2_true)
hold on 
plot(t, theta_crls(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLS with mu = 0.95')

subplot(3,2,5)
plot(t,theta1_true)
hold on 
plot(t, theta_vff(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'VFF-RLS')

subplot(3,2,6)
plot(t,theta2_true)
hold on 
plot(t, theta_vff(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'VFF-RLS')

% Evaluation metrics
% --- MSE ---
mse_theta1_rls = mean((theta_rls(:,1) - theta1_true).^2);
mse_theta1_vff = mean((theta_vff(:,1) - theta1_true).^2);

mse_theta2_rls = mean((theta_rls(:,2) - theta2_true).^2);
mse_theta2_vff = mean((theta_vff(:,2) - theta2_true).^2);

% --- Estimation variance ---
var_theta1_rls = var(theta_rls(:,1));
var_theta1_vff = var(theta_vff(:,1));

var_theta2_rls = var(theta_rls(:,2));
var_theta2_vff = var(theta_vff(:,2));


% % --- Settling time after t = 400 ---
% t_start = 400;
% eps = 0.05;  % Error treshold
% 
% % RLS settling time
% theta1_rls = theta_rls(:,1);
% settle1_rls = find(abs(theta1_rls(t_start:end) - theta1_true(t_start:end)) < eps, 1);
% if ~isempty(settle1_rls)
%     settling_time_rls1 = settle1_rls; % steps after t=400
% else
%     settling_time_rls1 = NaN;
% end
% 
% theta2_rls = theta_rls(:,2);
% settle2_rls = find(abs(theta2_rls(t_start:end) - theta2_true(t_start:end)) < eps, 1);
% if ~isempty(settle2_rls)
%     settling_time_rls2 = settle2_rls; % steps after t=400
% else
%     settling_time_rls2 = NaN;
% end
% 
% % VFF-RLS settling time
% theta1_vff = theta_vff(:,1);
% settle1_vff = find(abs(theta1_vff(t_start:end) - theta1_true(t_start:end)) < eps, 1);
% if ~isempty(settle1_vff)
%     settling_time_vff1 = settle1_vff;
% else
%     settling_time_vff1 = NaN;
% end
% 
% theta2_vff = theta_vff(:,2);
% settle2_vff = find(abs(theta2_vff(t_start:end) - theta2_true(t_start:end)) < eps, 1);
% if ~isempty(settle2_vff)
%     settling_time_vff2 = settle2_vff;
% else
%     settling_time_vff2 = NaN;
% end

% --- Misalignment ---
theta_true = [theta1_true, theta2_true];

misalign_rls = 10*log10(sum((theta_rls - theta_true).^2) ./ sum(theta_true.^2));
misalign_vff = 10*log10(sum((theta_vff - theta_true).^2) ./ sum(theta_true.^2));

avg_misalign_rls = mean(misalign_rls);
avg_misalign_vff = mean(misalign_vff);


% --- Print results ---
fprintf('--- Quantitative Comparison ---\n');
fprintf('Mean Square Error θ1 | RLS: %.4f    | VFF-RLS: %.4f\n', mse_theta1_rls, mse_theta1_vff);
fprintf('Mean Square Error θ2 | RLS: %.4f    | VFF-RLS: %.4f\n', mse_theta2_rls, mse_theta2_vff);
fprintf('Variance θ1          | RLS: %.4f    | VFF-RLS: %.4f\n', var_theta1_rls, var_theta1_vff);
fprintf('Variance θ2          | RLS: %.4f    | VFF-RLS: %.4f\n', var_theta2_rls, var_theta2_vff);
fprintf('Avg Misalignment(dB) | RLS: %.2f dB | VFF-RLS: %.2f dB\n', avg_misalign_rls, avg_misalign_vff);
% fprintf('Settling Time (θ1 after t=400)\n');
% fprintf('          RLS: %d steps\n', settling_time_rls1);
% fprintf('      VFF-RLS: %d steps\n', settling_time_vff1);
% fprintf('Settling Time (θ2)\n');
% fprintf('          RLS: %d steps\n', settling_time_rls2);
% fprintf('      VFF-RLS: %d steps\n', settling_time_vff2);


figure;
plot(misalign_rls, 'b'); hold on;
plot(misalign_vff, 'r');
legend('RLS', 'VFF-RLS');
xlabel('Time step'); ylabel('Misalignment (dB)');
title('Misalignment over Time');
grid on;


% Table creation
algorithms = {'RLS'; 'VFF-RLS'};
MSE_theta1 = [mse_theta1_rls; mse_theta1_vff];
MSE_theta2 = [mse_theta2_rls; mse_theta2_vff];
Variance_theta1 = [var_theta1_rls; var_theta1_vff];
Variance_theta2 = [var_theta2_rls; var_theta2_vff];
Misalignment_dB = [avg_misalign_rls; avg_misalign_vff];
% SettlingTime_theta1 = [settling_time_rls1; settling_time_vff1];
% SettlingTime_theta2 = [settling_time_rls2; settling_time_vff2];
T = table(MSE_theta1, MSE_theta2, Variance_theta1, Variance_theta2, Misalignment_dB, ...
          'RowNames', algorithms);
T = rows2vars(T);
writetable(T, 'comparison_results.csv', 'WriteRowNames', true);
