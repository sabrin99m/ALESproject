clearvars 
close all

rng('default');

% Simulation parameters
N = 45000;                             % Number of iterations
L = 64;                                % Filter length
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;

% Generate input signal (speech)
[x, Fs] = audioread('harvard.wav');
x = x(:,1);               % Usa solo primo canale
x = resample(x, 8000, Fs);% Porta a 8 kHz se necessario
x = x / max(abs(x));      % Normalizzazione

U = zeros(N, L);
for i = L:N
    U(i,:) = x(i:-1:i-L+1)';
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

% Generate noise signal(speech)
Fs = 8000;                     % Frequenza di campionamento 
start_idx = 22500;
end_idx = 33750;
snr_dB = 15;                   % SNR desiderato in dB
segment_len = end_idx - start_idx + 1;

% Segmento da corrompere
segment = d_noiseless(start_idx:end_idx);

% === 1. Genera segnale vocale sintetico (speech-like) ===
t = (0:segment_len-1)' / Fs;

% Formanti vocali (approssimazione per suono /a/)
f1 = 700; f2 = 1200; f3 = 2600;

speech_like = sin(2*pi*f1*t).*exp(-8*t) + ...
              0.5*sin(2*pi*f2*t).*exp(-6*t) + ...
              0.3*sin(2*pi*f3*t).*exp(-4*t);

% Normalizza
speech_like = speech_like / max(abs(speech_like));

% === 2. Scala per ottenere SNR = 15 dB ===
P_signal = mean(segment.^2);
P_noise_target = P_signal / 10^(snr_dB / 10);
P_noise_actual = mean(speech_like.^2);

scale_factor = sqrt(P_noise_target / P_noise_actual);
speech_noise_scaled = scale_factor * speech_like;

SNR_linear = 10^(SNR / 20);            % SNR linearization
var_d = var(d_noiseless);              % power of the desired noiseless signal
sigma_v = sqrt(var_d / SNR_linear);    % power of the system noise

v = zeros(N,1);
v(start_idx:end_idx) = speech_noise_scaled;
d = d_noiseless + v;

% Plain RLS (form 3) 
P = 10^6*eye(L);
theta_rls = zeros(N,L)';
lambda = 1-1/(10*L);
for i = L:N
    phi_i = U(i,:)';
    b = lambda + phi_i' * P * phi_i;
    K = (P * phi_i) / b;
    e = d(i) - phi_i' * theta_rls(:,i-1);
    theta_rls(:,i) = theta_rls(:,i-1) + K * e;
    P = (1/lambda)*(P - K * phi_i' * P);
end

theta_rls = theta_rls';


% GVFF
%___________________



% Variable Forgetting Factor RLS 

P = 10^6*eye(L);
theta_vff = zeros(N,L)';
lambda_max = 0.999999;                 % Maximum forgetting factor
gamma = 1.5;                           % lamba(n) = lambda_max when sigma_e(n) <= gamma*sigma_v(n)
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates

K_alpha = 6;
K_beta = 3*K_alpha;
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

% Evoluzione di lamnda nel tempo
figure
plot(L:N, lambda_n(L:N), 'r', 'DisplayName', 'VFF-RLS')
xlabel('Tempo [campioni]');
ylabel('lambda(n)');
legend;
title('Evoluzione di lambda(n)');
grid on;

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