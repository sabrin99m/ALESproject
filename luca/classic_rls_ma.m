function [h_hist] = classic_rls_ma(x,d,delta,lambda)
% x[N,1]
% P[L,L]
[N,kb] = size(x);
P = delta * eye(kb);
h_est = zeros(kb, 1); 
% Classic RLS
h_hist = zeros(N,kb);
y_est = zeros(N,1);
beta_hist = zeros(N,1);
P_hist = zeros(N,kb,kb);
for n = 1:N
    x_n = x(n,:)';
    beta = lambda + x_n'*P*x_n;
    P = (1/lambda)*(P - (1/beta) * P * (x_n * x_n') * P);
    k_n = P*x_n;
    y_est(n) = x_n'*h_est;
    e_n = d(n) - x_n'*h_est;
    h_est = h_est + k_n*e_n;
    h_hist(n,:) = h_est.';
    beta_hist(n,:) = beta;
    P_hist(n,:,:) = P;
end

figure
plot(beta_hist)
figure
plot(P_hist(:,1,1))

end

