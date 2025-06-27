function plotParameter(idx,real,est_classic,est_vff,est_gvff)
figure;
plot(real(:,idx), 'Color','black', 'DisplayName', 'Theta '+string(idx)+' Real');
hold on;
plot(est_classic(:,idx), 'b', 'DisplayName', 'Theta '+string(idx)+' Classic RLS');
plot(est_vff(:,idx), 'r', 'DisplayName', 'Theta '+string(idx)+' VFF')
plot(est_gvff(:,idx), 'g', 'DisplayName', 'Theta '+string(idx)+' GVFF')
legend;
title("RLS h"+string(idx)+" estimations")
grid on
ylim([-2,2])
end