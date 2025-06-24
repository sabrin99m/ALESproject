function [noise] = generate_noise(N,SNR_db,noiseless_signal)
% Generate noise signal
SNR_linear = 10^(SNR_db / 10);            % SNR linearization
var_d = var(noiseless_signal);              % power of the desired noiseless signal
sigma_v = sqrt(var_d / SNR_linear);    % power of the system noise

noise = sigma_v*randn(N,1);  
end

