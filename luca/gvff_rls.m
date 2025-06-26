function [w_hist, lambda_hist] = gvff_rls(u, y, delta, K_alpha, K_beta)
% x: vettore colonna (N x 1), segnale di input
% L: ordine del filtro
% d: vettore colonna (N x 1), segnale desiderato
% delta: inizializzazione matrice P
% mu: passo del gradiente
% lambda0: valore iniziale forgetting factor
% lambda_min, lambda_max: limiti di lambda

[N,L] = size(u);
theta = zeros(L, 1);              % pesi iniziali
P = delta * eye(L);           % matrice di correlazione inversa
mu = 0.1;
lambda_max = 0.999999;
lambda = lambda_max;
rho = 0;
rho_tilde = 0;
sigma_e = 0;  % power of the a priori error signal
sigma_v = 0;
sigma2_e_partial = 0;

alpha = 1 - 1/(K_alpha*L);
beta = 1 - 1/(K_beta*L);

w_hist = zeros(N,L);
lambda_hist = zeros(N, 1);
lambda_min_hist = zeros(N, 1);

for n = 1:N
    phi_n = u(n,:)';          % vettore colonna (L x 1)
    dn = y(n);                % segnale desiderato

        % Errore
    error = dn - theta' * phi_n;
    
    % Grad

    D = (L - 2)^2 * rho^2  - 8 * (L+2) * ((L+1) * rho_tilde + rho^2);
    if(D>0)
        lambda_min = ((L - 2) * rho + sqrt(D)) / (4 * ((L+1) * rho_tilde + rho^2));
    else
        lambda_min = 0;
    end

    rho = 1 + lambda*rho;
    rho_partial = rho;
    rho_tilde = 1 + lambda^2*rho_tilde;
    rho_tilde_partial = 2*lambda*rho_tilde;





    sigma2_e = alpha*sigma_e^2 + (1-alpha) * error^2;
    sigma_e = sqrt(sigma2_e);
    sigma2_v = beta*sigma_v^2 + (1-beta) * error^2;
    sigma_v = sqrt(sigma2_v);

    zeta = 1 - ((2*((L + 1)*rho_tilde + rho^2) - (L + 2)*rho)/(rho*((L+1)*rho_tilde+rho^2)));
    temp1 = 2/rho^2*rho_partial;
    temp2 = ((L + 1)*rho_tilde + rho^2)^2;
    temp3 = ((L + 1)* rho_tilde_partial + 2*rho*rho_partial);

    zeta_partial = temp1 - ((N + 2) /temp2 * temp3);
    h_partial = temp1 - (2 /temp2 * temp3);

    sigma2_e_partial = zeta*sigma2_e_partial + zeta_partial*sigma2_e + h_partial*sigma2_v;

    % Aggiorna lambda
    lambda = lambda - (mu * sigma2_e_partial / (1 - lambda));
    lambda = max(2*lambda_min, min(lambda_max, lambda));


    % Guadagno
    k = (P * phi_n) / (lambda + phi_n' * P * phi_n);

    % Aggiorna pesi
    theta = theta + k * error;

    % Aggiorna matrice P
    P = (P - k * phi_n' * P) / lambda;

    % Salva
    lambda_min_hist(n) = theta' * phi_n;
    w_hist(n,:) = theta;
    lambda_hist(n) = lambda;
end
plot(lambda_min_hist)
end