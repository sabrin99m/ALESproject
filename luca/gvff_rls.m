function [w_hist, lambda_hist] = gvff_rls(u, y, delta)
% x: vettore colonna (N x 1), segnale di input
% L: ordine del filtro
% d: vettore colonna (N x 1), segnale desiderato
% delta: inizializzazione matrice P
% mu: passo del gradiente
% lambda0: valore iniziale forgetting factor
% lambda_min, lambda_max: limiti di lambda

[N,L] = size(u);
w = zeros(L, 1);              % pesi iniziali
P = delta * eye(L);           % matrice di correlazione inversa
z = zeros(L, 1); 
mu = 0.1;
lambda_min = 0.5;
lambda0 = 0.99;
lambda_max = 0.999999;
lambda = lambda0;

w_hist = zeros(N,L);
lambda_hist = zeros(N, 1);

for n = 1:N
    xn = u(n,:)';          % vettore colonna (L x 1)
    dn = y(n);                % segnale desiderato

    % Errore
    error = dn - w' * xn;

    % Guadagno
    k = (P * xn) / (lambda + xn' * P * xn);

    % Aggiorna pesi
    w = w + k * error;

    % Aggiorna matrice P
    P = (P - k * xn' * P) / lambda;

    % Derivata ricorsiva z
    z = z + k * (error - xn' * z);

    % Gradiente
    gradJ = - error * (xn' * z);   % scalare

    % Aggiorna lambda
    lambda = lambda - mu * gradJ;
    lambda = max(lambda_min, min(lambda_max, lambda));

    % Salva
    w_hist(n,:) = w;
    lambda_hist(n) = lambda;
end
end