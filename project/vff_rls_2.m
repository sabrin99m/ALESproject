clc 
clearvars 
close all

rng('default');

% Simulation parameters
N = 45000;                             % Number of iterations
L = 64;                                % Filter length
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;


% Generate input signal (gaussian white noise)
sigma_u = 100;
u = sigma_u*randn(N,1);

U = zeros(N, L);
for i = L:N
    U(i,:) = u(i:-1:i-L+1)';
end

% Generate m1(k) system (first impulse response from ITU-T G168 Reccomendation) 
h_true = [...
    -436,  -829, -2797, -4208, -17968, -11215, ...
    46150, 34480, -10427,  9049, -1309, -6320, ...
    390, -8191, -1751, -6051, -3796, -4055, ...
    -3948, -2557, -3372, -1808, -2259, -1300, ...
    -1098,  -618,  -340,   -61,   323,   419, ...
    745,   716,   946,   880,  1014,   976, ...
    1033,  1091,  1053,  1042,   794,   831, ...
    899,   716,   390,   313,   304,   304, ...
    73,  -119,  -109,  -176,  -359,  -407, ...
    -512,  -580,  -704,  -618,  -685,  -791, ...
    -772,  -820,  -839,  -724,  ...
];

H = zeros(N,L);
for n = 1:N
    if n < 15000
        h_shift = h_true;
    elseif n < 30000
        h_shift = [zeros(1,4), h_true];
    else
        h_shift = [zeros(1,8), h_true];
    end
    H(n,:) = h_shift(1:L);  % salva h_shift a tempo n
end


d_noiseless = zeros(N,1);
for n = L:N
    d_noiseless(n) = U(n,:) * H(n,:)';
end

% Generate noise signal
SNR_linear = 10^(SNR / 20);            % SNR linearization
var_d = var(d_noiseless);              % power of the desired noiseless signal
sigma_v = sqrt(var_d / SNR_linear);    % power of the system noise

v = sigma_v*randn(N,1);                % noise signal 
d = d_noiseless + v;                   % desired signal

% Plain RLS (form 3) 
P = 10^6*eye(L);
theta_rls = zeros(N,L)';
lambda = 1-1/(3*L);
for i = L:N
    phi_i = U(i,:)';
    b = lambda + phi_i' * P * phi_i;
    K = (P * phi_i) / b;
    e = d(i) - phi_i' * theta_rls(:,i-1);
    theta_rls(:,i) = theta_rls(:,i-1) + K * e;
    P = (1/lambda)*(P - K * phi_i' * P);
end

theta_rls = theta_rls';

y_hat1 = zeros(N,1);

for i = L:N
    y_hat1(i) = U(i,:) * theta_rls(i-1,:)';
end

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
q = U(2,:) * P *U(2,:)';               
sigma_q = sqrt(var(q));                % prediction error associated with the input u(n)
sigma_e = sqrt(mean(d_noiseless.^2));  % power of the a priori error signal

lambda_n = ones(N, 1) * lambda_max;    % Variable Forgetting Factor initialization

for i = L:N
    phi_i = U(i,:)';
    e = d(i) - phi_i' * theta_vff(:,i-1);
    q_n = phi_i' * P * phi_i;
    K = (P * phi_i) / (lambda_n(i-1) + q_n);
    theta_vff(:,i) = theta_vff(:,i-1) + K * e;
    P = (1/lambda_n(i)) * (P - K * phi_i' * P);

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

y_hat2 = zeros(N,1);

for i = L:N
    y_hat2(i) = U(i,:) * theta_vff(i-1,:)';
end

% Plot results
figure
title('Plain RLS vs. VFF-RLS')
t = L:N;    % temporal line for plot

subplot(1,2,1)
plot(t,d(L:N))
hold on 
plot(t, y_hat1(L:N))
xlabel('time stamp')
ylabel('output')
legend('true system', 'estimated system (RLS)')

subplot(1,2,2)
plot(t,d(L:N))
hold on 
plot(t, y_hat2(L:N))
xlabel('time stamp')
ylabel('output')
legend('true system', 'estimated system (VFF)')

% Prealloca vettori per misalignment
misalignment_rls = zeros(N,1);
misalignment_vff = zeros(N,1);

for n = L:N
    % Aggiorna il vero h(n) in base alla shift logic
    if n < 15000
        h_shift = h_true;  % no shift
    elseif n < 30000
        h_shift = [zeros(1,4), h_true];  % shift di 4
    else
        h_shift = [zeros(1,8), h_true];  % shift di 8
    end
    h_shift = h_shift(1:L);  % prendi solo i primi L coefficienti

    % Calcola misalignment RLS
    err_rls = norm(theta_rls(n,:) - h_shift)^2;
    misalignment_rls(n) = 10 * log10(err_rls / norm(h_shift)^2);

    % Calcola misalignment VFF-RLS
    err_vff = norm(theta_vff(n,:) - h_shift)^2;
    misalignment_vff(n) = 10 * log10(err_vff / norm(h_shift)^2);
end

% Plotta
figure;
plot(L:N, misalignment_rls(L:N), 'b', 'DisplayName', 'Plain RLS');
hold on;
plot(L:N, misalignment_vff(L:N), 'r', 'DisplayName', 'VFF-RLS');
xlabel('Tempo [campioni]');
ylabel('Misalignment [dB]');
legend;
title('Misalignment vs tempo (in dB)');
grid on;

% --- Evaluation Metrics ---
% Mean Square Error (MSE) tra stima θ e vero H
mse_rls = mean(mean((theta_rls - H).^2));
mse_vff = mean(mean((theta_vff - H).^2));

% Varianza delle stime θ (calcolata su ogni coefficiente, poi mediata)
var_rls = mean(var(theta_rls));
var_vff = mean(var(theta_vff));

% Misalignment medio (media temporale su L:N)
avg_misalign_rls = mean(misalignment_rls(L:N));
avg_misalign_vff = mean(misalignment_vff(L:N));

% --- Print results ---
fprintf('\n--- Quantitative Comparison ---\n');
fprintf('Mean Square Error θ      | RLS: %.4f     | VFF-RLS: %.4f\n', mse_rls, mse_vff);
fprintf('Estimation Variance θ    | RLS: %.4f     | VFF-RLS: %.4f\n', var_rls, var_vff);
fprintf('Avg Misalignment (dB)    | RLS: %.2f dB  | VFF-RLS: %.2f dB\n', avg_misalign_rls, avg_misalign_vff);


% Table creation
% algorithms = {'RLS'; 'VFF-RLS'};
% MSE_theta1 = [mse_theta1_rls; mse_theta1_vff];
% MSE_theta2 = [mse_theta2_rls; mse_theta2_vff];
% Variance_theta1 = [var_theta1_rls; var_theta1_vff];
% Variance_theta2 = [var_theta2_rls; var_theta2_vff];
% Misalignment_dB = [avg_misalign_rls; avg_misalign_vff];
% % SettlingTime_theta1 = [settling_time_rls1; settling_time_vff1];
% % SettlingTime_theta2 = [settling_time_rls2; settling_time_vff2];
% T = table(MSE_theta1, MSE_theta2, Variance_theta1, Variance_theta2, Misalignment_dB, ...
%           'RowNames', algorithms);
% T = rows2vars(T);
% writetable(T, 'comparison_results.csv', 'WriteRowNames', true);
