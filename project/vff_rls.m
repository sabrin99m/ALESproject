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
% Adaptive RLS (lambda = 0.95)
% GVFF RLS
% Variable Forgetting Factor RLS 

% Plain RLS (form 3) 
P = 10^6*eye(L);
theta = zeros(N,L)';

for i = L:N
    b = 1 + u(i,:)*P*u(i,:)';
    P = P - (1/b)*P*u(i,:)'*u(i,:)*P;
    K = P*u(i,:)';
    e = y(i) - u(i,:)*theta(:,i-1);
    theta(:,i) = theta(:,i-1) + K*e;
end

theta = theta';

% Adaptive RLS (lambda = 0.95)

P = 10^6*eye(L);
theta2 = zeros(N,L)';
lambda = .95;
for i = L:N
    b = lambda + u(i,:)*P*u(i,:)';
    P = (1/lambda)*(P - (1/b)*P*u(i,:)'*u(i,:)*P);
    K = P*u(i,:)';
    e = y(i) - u(i,:)*theta2(:,i-1);
    theta2(:,i) = theta2(:,i-1) + K*e;
end

theta2 = theta2';


% GVFF
%___________________



% Variable Forgetting Factor RLS 

P = 10^6*eye(L);
theta3 = zeros(N,L)';
lambda_max = 0.999999;                 % Maximum forgetting factor
gamma = 1.5;                           % lamba(n) = lambda_max when sigma_e(n) <= gamma*sigma_v(n)
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates

K_alpha = 2;
K_beta = 5*K_alpha;
alpha = 1 - 1/(K_alpha*L);
beta = 1 - 1/(K_beta*L);

sigma_q = sqrt(L);                     % stima iniziale
sigma_e = sqrt(mean(y_noiseless.^2));  % power of the a priori error signal

lambda_n = ones(N, 1) * lambda_max;    % Variable Forgetting Factor initialization

for i = 2:N
    e = y(i) - u(i,:)*theta3(:,i-1);                            % A priori error
    K = (P * u(i,:)') / (lambda_n(i-1) + u(i,:) * P * u(i,:)');  % Kalman gain vector
    theta3(:,i) = theta3(:,i-1) + K*e;
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

theta3 = theta3';


% Plot results
figure
t = 1:N;    % temporal line for plot

subplot(3,2,1)
plot(t,theta1_true)
hold on 
plot(t, theta(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLSplain')

subplot(3,2,2)
plot(t,theta2_true)
hold on 
plot(t, theta(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLSplain')

subplot(3,2,3)
plot(t,theta1_true)
hold on 
plot(t, theta2(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLS with mu = 0.95')

subplot(3,2,4)
plot(t,theta2_true)
hold on 
plot(t, theta2(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'RLS with mu = 0.95')

subplot(3,2,5)
plot(t,theta1_true)
hold on 
plot(t, theta3(:,1))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'VFF-RLS')

subplot(3,2,6)
plot(t,theta2_true)
hold on 
plot(t, theta3(:,2))
xlabel('time stamp')
ylabel('theta')
legend('theta-true', 'VFF-RLS')
