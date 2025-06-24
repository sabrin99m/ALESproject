function plotParameter(idx,real,est)
figure;
plot(est(:,idx), 'b', 'DisplayName', 'h'+string(idx)+'_est');
hold on;
plot(real(:,idx), 'g', 'DisplayName', 'h'+string(idx)+'_real','LineStyle','--');
legend('Estimated', 'Real', 'Location', 'best');
title("Normal RLS h"+string(idx)+" estimations")
%ylim([-2,2])
end