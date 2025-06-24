function [misalignment] = puntual_misalignment(estimated,real)
misalignment = 20*log10(norm(real - estimated) / norm(real));
end

