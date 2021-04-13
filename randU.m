function U = randU(n)
    X = (randn(n) + 1j*randn(n))/sqrt(2);
    [Q,R] = qr(X);
    R = diag(diag(R) ./ abs(diag(R)));
    U = Q*R;
end