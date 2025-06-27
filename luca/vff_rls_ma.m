function [h_hist,lambda_hist] = vff_rls_ma(u,y,y_noiseless,delta,K_alpha,K_beta)
% RLS Initialization
[N,L] = size(u);
h_est = zeros(L, 1);                   % Adaptive filter
gamma = 1.5;   
P = eye(L) * delta;        % Inverse of the input auto-correlation matrix
xi = 10^-8;                            % Small constant to avoid division per zero in forgetting factor updates
lambda_max = 0.999999;  
lambda = lambda_max;    % Variable Forgetting Factor 




alpha = 1 - 1/(K_alpha*L);
beta = 1 - 1/(K_beta*L);

q = u(2,:) * P *u(2,:)';               
sigma_q = sqrt(var(q));                % prediction error associated with the input u(n)
sigma_e = sqrt(mean(y_noiseless.^2));  % power of the a priori error signal
sigma_v = 1;

% Storage for tracking
h_hist = zeros(N, L);
lambda_hist = zeros(N,1);
% RLS with Variable Forgetting Factor
for n = 1:N
    x_n = u(n,:)';  % Input vector

    e_n = y(n) - h_est' * x_n; % A priori error
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
    
    if(sigma_e > gamma*sigma_v)
        lambda = min((sigma_q*sigma_v)/(xi + abs(sigma_e - sigma_v)), lambda_max);
    else
        lambda = lambda_max;
    end
    
    % Store metrics
    h_hist(n,:) = h_est;
    lambda_hist(n) = lambda;
end

end