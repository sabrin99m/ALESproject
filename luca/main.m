clc 
clearvars 
close all

rng('default');
set(0, 'DefaultFigureVisible', 'on');

%% Initialization

% Simulation parameters
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;
P_delta = 10^6;                        % Initial value of P (Suggested high value)

%% Course system realizations
N = 1000;
L=2;
e = rand(N,1)*sqrt(1);
h = course_theta(N);
x = course_input(N);
y = ma_realizations(h,x,e);
%% RLS
h_classic = normal_rls_procedure(x,y,h,P_delta);
h_vff = vff_rls_procedure(x,y,h,P_delta);
h_gvff = gvff_rls_procedure(x,y,h,P_delta);
%% Paper system realization
N = 45000;
L=64;
h = paper_theta(N,L);
x = paper_input_w(N,L);
v = ma_realizations(h,x,zeros(N,1));
e = generate_noise(N,SNR,v);
y = ma_realizations(h,x,e);
%% RLS
h_classic = normal_rls_procedure(x,y,h,P_delta);
h_vff = vff_rls_procedure(x,y,h,P_delta);
h_gvff = gvff_rls_procedure(x,y,h,P_delta);

%% Plots
plotParameter(1,h,h_classic)
plotParameter(1,h,h_vff)
plotParameter(1,h,h_gvff)
plotParameter(2,h,h_classic)
plotParameter(2,h,h_vff)
plotParameter(2,h,h_gvff)
%% Helper functions
function h_rls = normal_rls_procedure(x,y,h,P_delta)
h_rls = classic_rls_ma(x,y,P_delta,0.95);
misalignment(h_rls,h);
end

function [h_rls,lambda] = vff_rls_procedure(x,y,h,P_delta)
[h_rls,lambda] = vff_rls_ma(x,y,P_delta);
misalignment(h_rls,h);
figure;
plot(lambda)
end

function [h_rls,lambda] = gvff_rls_procedure(x,y,h,P_delta)
[h_rls,lambda] = gvff_rls(x,y,P_delta);
misalignment(h_rls,h);
figure;
plot(lambda)
end

function x = course_input(N)
sigma_u1 = sqrt(100); % Deviazione standard per u1
sigma_u2 = sqrt(100); % Deviazione standard per u2

w1 = rand(N,1)*sigma_u1;
w2 = rand(N,1)*sigma_u2;
x=[w1,w2];
end

function x = paper_input_w(N,L)
sigma_w = 100;
w = sigma_w*randn(N,1);
x = get_lags(w,L);
end

function x = paper_input_ar1(N,L)
% Filtro AR(1): H(z) = 1 / (1 - a*z^-1)
w = rand(N,1)*sqrt(100);
a = 0.9;                % Coeff
b = 1;                  % Numeratore
a_vec = [1, -a];        % Denominatore: 1 - 0.9*z^-1

x_ar1 = filter(b, a_vec, w);
x = get_lags(x_ar1,L);
end

function h = course_theta(N)
theta_2 = ones(N,1);
theta_1 = ones(N,1);
for n = 400:N
    theta_1(n)=sin(0.01*pi*n + pi/2);       
end
h = [theta_1,theta_2];
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

function y = ma_realizations(theta, phi, epsilon)
N = size(epsilon,1);
y = zeros(N,1);
for n = 1:N
    y(n) = phi(n,:)*theta(n,:)' + epsilon(n);
end

end