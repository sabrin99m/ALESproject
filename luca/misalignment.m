function misalignment = misalignment(est,real)
arguments
    est (:,:) double
    real (:,:) double
end
N = size(est,1);
misalignment = zeros(N,1);
for n = 1:N
    misalignment(n) = puntual_misalignment(est(n,:),real(n,:));
end
end