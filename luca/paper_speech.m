%% INIT
clc 
clearvars 
close all

rng('default');

%% PARAMETERS
N = 45000;
L=64;
SNR = ones(N,1)*20;                              % dB, snr = sigma_d^2/sigma_v^2;
P_delta = 10^6; 
lambda_rls = 1-(1/(10*L));

%% SYSTEM
h = paper_theta(N,L);
x = paper_input(N,L);
v = ma_realizations(h,x,zeros(N,1));
e = generate_noise(N,SNR,v);
y = ma_realizations(h,x,e);
d = ma_realizations(h,x,zeros(N,1));

%% RLS
[h_classic,mis_classic] = normal_rls_procedure(x,d,h,P_delta,lambda_rls);
[h_vff, lambda_vff, mis_vff] = vff_rls_procedure(x,d,d,h,P_delta,6,18);
[h_gvff, lambda_gvff, mis_gvff] = gvff_rls_procedure(x,d,h,P_delta);

%% PLOTS
plot_mis_and_lambda(mis_classic,mis_vff,mis_gvff,lambda_vff,lambda_gvff)

function x = paper_input(N,L)
[x, Fs] = audioread('harvard.wav');
x = x(:,1);               % Usa solo primo canale
x = resample(x, 8000, Fs);% Porta a 8 kHz se necessario
x = x / max(abs(x));      % Normalizzazione
x = get_lags(x(1:N),L);
end

function h = paper_theta(N,L)
h_init = [...
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

h = zeros(N,L);
h_shift = [zeros(1,4), h_init(1:end-4)];
h_shift_shift = [zeros(1,4), h_shift(1:end-4)];

for n = 1:N
    if n < 15000
        h(n,:) = h_init;  % no shift
    elseif n < 30000
        h(n,:) = h_shift;  % shift of 4
    else
        h(n,:) = h_shift_shift;  % shift of 4 again
    end
end
end