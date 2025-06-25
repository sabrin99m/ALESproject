function [h_rls,mis] = normal_rls_procedure(x,y,h,P_delta, lambda)
h_rls = classic_rls_ma(x,y,P_delta,lambda);
mis = misalignment(h_rls,h);
end

