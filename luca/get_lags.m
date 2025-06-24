function X_lag=get_lags(x,L)
N = length(x);
    X_lag = zeros(N, L);  % inizializza la matrice con zeri

    for i = 1:N
        idx_start = max(1, i - L + 1);  % indice di inizio, senza scendere sotto 1
        valori = x(idx_start:i);       % prendi i valori disponibili
        X_lag(i, end-length(valori)+1:end) = valori;  % riempi da destra
    end
end