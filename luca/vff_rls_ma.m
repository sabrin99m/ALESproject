function [h_hist,lambda_hist] = vff_rls_ma(x,d,delta)
% RLS Initialization
[N,kb] = size(x);
h_est = zeros(kb, 1);                   % Adaptive filter
P = eye(kb) * delta;        % Inverse of the input auto-correlation matrix
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates
lambda_max = 0.999999;  
lambda = lambda_max;    % Variable Forgetting Factor 

K_alpha = 2;
K_beta = 5*K_alpha;
alpha = 1 - 1/(K_alpha*kb);
beta = 1 - 1/(K_beta*kb);

sigma_q = 0;                     % stima iniziale
sigma_e = 0;                 % power of the a priori error signal
sigma_v = 0;

% Storage for tracking
h_hist = zeros(N, kb);
lambda_hist = zeros(N,1);
% RLS with Variable Forgetting Factor
for n = 1:N
    x_n = x(n,:)';  % Input vector
    e_n = d(n) - h_est' * x_n; % A priori error
    k_n = (P * x_n) / (lambda  + x_n' * P * x_n);  % Kalman gain vector
    h_est = h_est + k_n * e_n;
    q_n = x_n' * P *x_n;
    P = (P - (k_n * x_n' * P)) / lambda;  %inverse of input correlation matrix
    
    % Update forgetting factor
    sigma2_e = alpha*sigma_e^2 + (1-alpha)*e_n^2;
    sigma_e = sqrt(sigma2_e);
    sigma2_q = alpha*sigma_q^2 + (1-alpha)*q_n^2;
    sigma_q = sqrt(sigma2_q);
    sigma2_v = beta*sigma_v^2 + (1-beta)*e_n^2;
    sigma_v = sqrt(sigma2_v);

    lambda = min((sigma_q*sigma_v)/(xi + abs(sigma_e - sigma_v)), lambda_max);
    
    % Store metrics
    h_hist(n,:) = h_est;
    lambda_hist(n) = lambda;
end

end