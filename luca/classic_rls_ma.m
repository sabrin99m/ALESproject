function [theta_hist] = classic_rls_ma(u,y,delta,lambda)
% x[N,1]
% P[L,L]
[N,kb] = size(u);
P = delta * eye(kb);
theta_est = zeros(kb, 1); 
% Classic RLS
theta_hist = zeros(N,kb);
y_est = zeros(N,1);
for n = 1:N
    phi = u(n,:)';
    beta = lambda + phi'*P*phi;
    P = (1/lambda)*(P - (1/beta) * P * (phi * phi') * P);
    k = P*phi;
    y_est(n) = phi'*theta_est;
    error = y(n) - phi'*theta_est;
    theta_est = theta_est + k*error;
    

    theta_hist(n,:) = theta_est.';
end
end

