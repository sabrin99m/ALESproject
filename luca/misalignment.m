function misalignment = misalignment(est,real, doPlot)
arguments
    est (:,:) double
    real (:,:) double
    doPlot (1,1) logical = true
end
N = size(est,1);
misalignment = zeros(N,1);
for n = 1:N
    p_mis = puntual_misalignment(est(n,:),real(n,:));
    misalignment(n) = puntual_misalignment(est(n,:),real(n,:));
end
if doPlot
figure;
plot(misalignment);
title("Misalignment");
end
end