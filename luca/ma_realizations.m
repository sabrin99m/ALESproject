function y = ma_realizations(theta, phi, epsilon)
N = size(epsilon,1);
y = zeros(N,1);
for n = 1:N
    y(n) = phi(n,:)*theta(n,:)' + epsilon(n);
end
end