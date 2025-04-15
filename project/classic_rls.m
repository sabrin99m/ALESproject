clc 
clearvars 
close all

rng('default');

% Simulation parameters
N = 45000;                             % Number of iterations
L = 64;                                % Filter length
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;


% Generate input signal (gaussian white noise)
sigma_x = 100;
x = sigma_x*randn(N,1);


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
        h_shift = [zeros(1,4), h_true];  % shift of 4
    end

    h_shift = h_shift(1:L);
end


% Desired noiseless signal
d_noiseless = filter(h_shift,1,x);

% Generate noise signal
SNR_linear = 10^(SNR / 20);            % SNR linearization
var_d = var(d_noiseless);              % power of the desired noiseless signal
sigma_v = sqrt(var_d / SNR_linear);    % power of the system noise

v = sigma_v*randn(N,1);                % noise signal 
d = d_noiseless + v;                   % desired signal

% RLS Initialization
h_est = zeros(L, 1);                   % Adaptive filter
P = eye(L) * 10^6;                     % Inverse of the input auto-correlation matrix
lambda = 1 - 1/(3*L);                  % Constant Forgetting Factor 


% Storage for tracking
misalignment = zeros(N, 1);

% RLS with Variable Forgetting Factor
for n = L:N
    x_n = x(n-L+1:n)';  % Input vector
    beta = lambda + x_n*P*x_n';
    P = (1/lambda)*(P - (1/beta)*P*x_n'*x_n*P);
    k_n = P*x_n';
    e_n = d(n) - x_n*h_est;
    h_est = h_est + k_n*e_n;
    
    % Store metrics
    misalignment(n) = 20*log10(norm(h_shift - h_est) / norm(h_shift));
end

% Plot results
figure;
subplot(2,1,1);
plot(misalignment(65:N), 'LineWidth', 1.5);
xlabel('Iterations');
ylabel('Misalignment [dB]');
xlim([0,N])
ylim([-50,30])
title('Misalignment over Time');
grid on;

subplot(2,1,2);
plot(ones(N, 1)*lambda, 'LineWidth', 1.5);
xlabel('Iterations');
ylabel('Forgetting Factor \lambda(n)');
title('Adaptive Forgetting Factor');
grid on;
