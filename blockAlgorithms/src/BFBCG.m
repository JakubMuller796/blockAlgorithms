function [x,bf_omega,numrank,conv] = BFBCG(A,b,x_0,type,tol,rankrevtol,maxit,xex)

tic
bf_omega = [];
T = [];
[n,m] = size(b);
conv = [0,0];

x = x_0;

r = b - A*x_0;
if m > 1
    [p,~,~] = qr(r,"econ");
else
    p = r/norm(r);
end

numrank = [size(p,2),0];
curr_size = size(p,2);

den = norm(b,'fro');

for i = 1:maxit
    Q = A*p;

    alpha = (p'*Q)\(p'*r);
    x = x + p*alpha;

    r = r - Q*alpha;
    
    bf_omega = [bf_omega trace((x-xex)'*A*(x-xex))];
    
    if tol > 0 && conv(1) == 0
        res = norm(r,'fro') / den;
        if res < tol
            % fprintf("Convergence in %d\n",i)
            conv(1) = toc;
            conv(2) = i;
        end
    end

    beta = -(p'*Q)\(Q'*r);
    
    if m > 1
        if type == "svd"
            %%% SVD
            [p,G,~] = svd(r + p*beta,"econ");
            d = abs(diag(G));
            rank = sum(d > rankrevtol*d(1));
            if rank ~= curr_size;
                numrank = [numrank;[rank, i]];
                curr_size = rank;
            end
            p = p(:,1:rank);
        elseif type == "qr"
            %%% QR
            [p,G,~] = qr(r + p*beta,"econ","vector");
            d = abs(diag(G(1:size(p,2),1:size(p,2))));
            rank = sum(d > rankrevtol*d(1));
            if rank ~= curr_size;
                numrank = [numrank;[rank, i]];
                curr_size = rank;
            end
            p = p(:,1:rank);
        end
    else
        %%% Vector
        p = (r + p*beta)/norm(r + p*beta);
    end
end