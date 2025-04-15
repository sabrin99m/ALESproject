clc 
clearvars 
close all

rng('default');

% Simulation parameters
N = 45000;                             % Number of iterations
L = 64;                                % Filter length
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates
lambda_max = 0.999999;                 % Maximum forgetting factor
gamma = 1.5;                           % lamba(n) = lambda_max when sigma_e(n) <= gamma*sigma_v(n)
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;

K_alpha = 2;
K_beta = 5*K_alpha;
alpha = 1 - 1/(K_alpha*L);
beta = 1 - 1/(K_beta*L);

% Generate input signal (gaussian white noise)
sigma_x = 100;
x = sigma_x*randn(N,1);
x_n = x(N:-1:N-L+1,1);

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

for n = 1:N
    if n < 15000
    h_shift = [h_true];  % no shift
    elseif n < 30000
        h_shift = [zeros(1,4), h_true];  % shift of 4
    else
        h_shift = [zeros(1,8), h_true];  % shift of 8
    end

    h_shift = h_shift(1:L);
end

% Desired noiseless signal
d_noiseless = filter(h_shift, 1, x);

% Generate noise signal
SNR_linear = 10^(SNR / 20);            % SNR linearization
var_d = var(d_noiseless);              % power of the desired noiseless signal
sigma_v = sqrt(var_d / SNR_linear);    % power of the system noise

v = sigma_v*randn(N,1);                % noise signal 
d = d_noiseless + v;                   % desired signal

% RLS Initialization
h_est = zeros(L, 1);                   % Adaptive filter
P = eye(L) * 10e6;                     % Inverse of the input auto-correlation matrix
lambda_n = ones(N, 1) * lambda_max;    % Variable Forgetting Factor 

sigma_q = sqrt(L);                     % stima iniziale
sigma_e = sqrt(var_d);                 % power of the a priori error signal

% Storage for tracking
misalignment = zeros(N, 1);

% RLS with Variable Forgetting Factor
for n = L:N
    x_n = x(n-L+1:n);  % Input vector
    e_n = d(n) - h_est' * x_n; % A priori error
    k_n = (P * x_n) / (lambda_n(n-1) + x_n' * P * x_n);  % Kalman gain vector
    h_est = h_est + k_n * e_n;
    q_n = x_n' * P *x_n;
    P = (P - (k_n * x_n' * P)) / lambda_n(n-1);  %inverse of input correlation matrix
    
    % Update forgetting factor
    sigma2_e = alpha*sigma_e^2 + (1-alpha)*e_n^2;
    sigma_e = sqrt(sigma2_e);
    sigma2_q = alpha*sigma_q^2 + (1-alpha)*q_n^2;
    sigma_q = sqrt(sigma2_q);
    sigma2_v = beta*sigma_v^2 + (1-beta)*e_n^2;
    sigma_v = sqrt(sigma2_v);

    lambda_n(n) = min((sigma_q*sigma_v)/(xi + abs(sigma_e - sigma_v)), lambda_max);
    
    % Store metrics
    misalignment(n) = 20*log10(norm(h_shift - h_est) / norm(h_shift));
end

% Plot results
figure;
subplot(2,1,1);
plot(misalignment(501:N), 'LineWidth', 1.5);
xlabel('Iterations');
ylabel('Misalignment [dB]');
xlim([0,N])
ylim([-50,30])
title('Misalignment over Time');
grid on;

% subplot(2,1,2);
% plot(lambda_n, 'LineWidth', 1.5);
% xlabel('Iterations');
% ylabel('Forgetting Factor \lambda(n)');
% title('Adaptive Forgetting Factor');
% grid on;
