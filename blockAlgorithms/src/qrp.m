function [Q,R] = qrp(A)

    [Q,R] = qr(A,"econ");
    I = speye(size(R));
    I = spdiags([sign(diag(R))], 0, I);
    Q = Q * I;
    R = I * R;
