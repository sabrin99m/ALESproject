%% INIT
clc 
clearvars 
close all

rng('default');

%% PARAMETERS
SNR = 20;                              % dB, snr = sigma_d^2/sigma_v^2;
P_delta = 10^6; 
N = 1000;
L=2;
lambda_rls = 0.95;

%% SYSTEM
e = rand(N,1)*sqrt(1);
h = course_theta(N);
x = course_input(N);
y = ma_realizations(h,x,e);
d = ma_realizations(h,x,zeros(N,1));

%% RLS
[h_classic,mis_classic] = normal_rls_procedure(x,y,h,P_delta,lambda_rls);
[h_vff, lambda_vff, mis_vff] = vff_rls_procedure(x,y,d,h,P_delta,2,10);
[h_gvff, lambda_gvff, mis_gvff] = gvff_rls_procedure(x,y,h,P_delta);

%% PLOTS
for i = 1:L
plotParameter(i,h,h_classic,h_vff,h_gvff)
end
plot_mis_and_lambda(mis_classic,mis_vff,mis_gvff,lambda_vff,lambda_gvff)

function x = course_input(N)
sigma_u1 = sqrt(100); % Deviazione standard per u1
sigma_u2 = sqrt(100); % Deviazione standard per u2

w1 = rand(N,1)*sigma_u1;
w2 = rand(N,1)*sigma_u2;
x=[w1,w2];
end

function h = course_theta(N)
theta_2 = ones(N,1);
theta_1 = ones(N,1);
for n = 400:N
    theta_1(n)=sin(0.01*pi*n + pi/2);       
end
h = [theta_1,theta_2];
end