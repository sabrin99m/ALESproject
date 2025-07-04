function [w_est, omega_hist] = gvff_rls(u, y, delta, K_alpha, K_beta, mu_omega, omega_max)
%GVFF_RLS_DYNAMIC Implementa l'algoritmo GVFF-RLS con limite inferiore di omega dinamico.
%
%   RIFERIMENTO:
%   S. H. Leung and C. F. So, "Gradient-Based Variable Forgetting Factor 
%   RLS Algorithm in Time-Varying Environments," IEEE Transactions on 
%   Signal Processing, vol. 53, no. 8, pp. 3141-3145, Aug. 2005.
%
%   SINTASSI:
%   [w_est, e, omega] = GVFF_RLS_dynamic(u, y, delta, K_alpha, K_beta, mu_omega, omega_max)
%
%   INPUT:
%   u            - Matrice degli ingressi [N_samples, L].
%   y            - Vettore colonna dell'uscita desiderata [N_samples, 1].
%   delta        - Valore di inizializzazione per la matrice P.
%   K_alpha      - Parametro alpha per la stima dell'MSE (Eq. 45).
%   K_beta       - Parametro beta per la stima della varianza del rumore (Eq. 46).
%   mu_omega     - Step-size (mu) per l'aggiornamento di omega (Eq. 44).
%   omega_max    - Limite superiore (ceiling) per il fattore di dimenticanza.
%
%   OUTPUT:
%   w_est        - Pesi del filtro stimati [L, N_samples].
%   e            - Errore a priori [N_samples, 1].
%   omega        - Evoluzione del fattore di dimenticanza [N_samples, 1].

%% --- Inizializzazione ---
[N, L] = size(u); % L è N nel paper
w = zeros(L, 1);
P = (1/delta) * eye(L);
mu_omega = 0.1;
omega_max = 0.999999;

% Inizializza omega al suo valore massimo, come da logica del paper
omega_val = omega_max;

% Variabili GVFF
sigma_e_sq = 1e-6;
sigma_eta_sq = 1e-6;
rho = 1;
rho_tilde = 1;
grad = 0;

% Pre-allocazione output
w_est = zeros(N,L);
e = zeros(N, 1);
omega_hist = zeros(N, 1);
eps = 1e-9; % Per stabilità numerica

%% --- Ciclo Principale ---
for n = 1:N
    x = u(n, :)';

    % --- Parte RLS Standard ---
    k = (P * x) / (omega_val + x' * P * x);
    e_n = y(n) - w' * x;
    w = w + k * e_n;
    P = (1/omega_val) * (P - k * x' * P);

    % --- Parte GVFF (calcolo del prossimo omega) ---
    sigma_e_sq_new = K_alpha * sigma_e_sq + (1-K_alpha) * e_n^2;
    sigma_eta_sq_new = K_beta * sigma_eta_sq + (1-K_beta) * e_n^2;
    rho_new = 1 + omega_val * rho;
    rho_tilde_new = 1 + omega_val^2 * rho_tilde;
    
    d_rho_d_omega = rho;
    d_rho_tilde_d_omega = 2 * omega_val * rho_tilde;
    denominator = (L+1) * rho_tilde_new + rho_new^2;
    common_term_d = ((L+1)*d_rho_tilde_d_omega + 2*rho_new*d_rho_d_omega);
    d_zeta_d_omega = (2/rho_new^2) * d_rho_d_omega - ((L+2)/denominator^2) * common_term_d;
    d_h_d_omega = -(2/rho_new^2) * d_rho_d_omega - (2/denominator^2) * common_term_d;
    zeta = 1 - 2/rho_new + (L+2)/denominator;
    
    grad_new = zeta * grad + d_zeta_d_omega * sigma_e_sq_new + d_h_d_omega * sigma_eta_sq_new;
    omega_next = omega_val - (mu_omega / (1 - omega_val + eps)) * grad_new;

    % --- Clipping con Limite Inferiore Dinamico (Eqs. 31, 32, 33) ---
    % 1. Calcolo del discriminante D. Usa rho e rho_tilde dello step precedente.
    D = (L-2)^2 * rho^2 - 8*(L+2)*((L+1)*rho_tilde + rho^2); % [cite: 197]
    
    % 2. Calcolo del limite inferiore dinamico omega_min = 2*omega#
    if D < 0
        % Se D < 0, la condizione di stabilità è sempre soddisfatta per omega > 0[cite: 196].
        % Il paper non definisce un floor in questo caso. Usiamo un valore
        % pratico e sicuro per evitare che omega scenda a valori troppo bassi.
        omega_min_dynamic = 0.5;
    else
        % Calcola omega# usando l'Eq. 33
        omega_hash_num = (L-2)*rho + sqrt(D);
        omega_hash_den = 4*((L+1)*rho_tilde + rho^2);
        omega_hash = omega_hash_num / (omega_hash_den + eps);
        omega_min_dynamic = 2 * omega_hash; % [cite: 257, 291]
    end
    
    % 3. Applica il clipping
    omega_next = max(omega_min_dynamic, min(omega_max, omega_next));
    
    % --- Aggiornamento e Salvataggio ---
    w_est(n, :) = w;
    e(n) = e_n;
    omega_hist(n) = omega_val;
    
    omega_val = omega_next;
    sigma_e_sq = sigma_e_sq_new;
    sigma_eta_sq = sigma_eta_sq_new;
    rho = rho_new;
    rho_tilde = rho_tilde_new;
    grad = grad_new;
end

end