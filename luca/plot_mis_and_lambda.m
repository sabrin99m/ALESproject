function plot_mis_and_lambda(mis_classic,mis_vff,mis_gvff,lambda_vff,lambda_gvff)
figure;
plot(mis_classic, 'b', 'DisplayName', 'Plain RLS');
hold on;
plot(mis_vff, 'r', 'DisplayName', 'VFF-RLS')
plot(mis_gvff, 'g', 'DisplayName', 'GVFF-RLS')
xlabel('Time [samples]');
ylabel('Misalignment [dB]');
legend;
title('Misalignment vs time (dB)');
grid on;
% ylim([-50,10])

figure;
plot(lambda_vff,'r', 'DisplayName', 'VFF-RLS');
hold on
xlabel('Time [samples]');
ylabel('lambda');
legend;
title('Evolution of lambda');
grid on;
plot(lambda_gvff,'b', 'DisplayName', 'GVFF-RLS');
end

