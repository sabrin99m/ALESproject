function [noise] = generate_noise(N,SNR_db,noiseless_signal)
% Generate noise signal
SNR_linear = zeros(N,1);
sigma_v = zeros(N,1);
var_d = var(noiseless_signal);
noise = randn(N,1);
for i = 1:N
SNR_linear(i) = 10^(SNR_db(i) / 10);
sigma_v(i) = sqrt(var_d / SNR_linear(i));
noise(i) = sigma_v(i)*noise(i);
end 
end

