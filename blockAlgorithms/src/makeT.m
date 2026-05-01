function T = makeT(T, alpha, beta)
[n,m] = size(alpha);

if isempty(T)
    if isempty(beta)
        T = alpha;
    else
        T = zeros(2*m);
        T(1:m,1:m) = alpha;
        T(1:m,m+1:2*m) = beta';
        T(m+1:2*m,1:m) = beta;
    end
    return;
end

n_old = size(T,1);

if isempty(beta)
    T(n_old-m+1:n_old, n_old-m+1:n_old) = alpha;
else
    Tnew = zeros(n_old + m);
    Tnew(1:n_old,1:n_old) = T;

    idx_curr = n_old - m + (1:m);
    idx_next = n_old + (1:m);

    % alpha
    Tnew(idx_curr, idx_curr) = alpha;

    % beta
    Tnew(idx_curr, idx_next) = beta';
    Tnew(idx_next, idx_curr) = beta;

    T = Tnew;
end
end