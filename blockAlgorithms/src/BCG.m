function [x,T,omega_bcg,numrank,conv] = BCG(A,b,x_0,maxit,tillConv,coef,reo,tol,rankrevtol,xex)

tic
conv = [0,0];
[n,m] = size(b);
omega_bcg = [];
numrank = [];
T = [];
den = norm(b,'fro');

x = x_0;
r = b - A*x_0;
p = r;

if reo, [V,~] = qrp(r); end

if coef
    [~,sigma] = qrp(r); 
    beta = sigma;
    theta = eye(size(r,2)); 
    ell = zeros(m,m);
end

numrank = [size(p,2),0];
curr_size = size(p,2);

for i = 1:maxit

    ga = (p'*A*p) \ (r'*r);
    x = x + p*ga;
    rnext = r - A*p*ga;

    if reo, [rnext,V] = reorth(rnext,V,reo); end 

    if isempty(xex) == 0
        omega_bcg = [omega_bcg trace((x-xex)'*A*(x-xex))];
    end

    if tol > 0 && conv(1) == 0
        res = norm(rnext,'fro') / den;
        if res < tol
            % fprintf("Convergence in %d\n",i)
            conv(1) = toc;
            conv(2) = i;
            if tillConv
                break
            end
        end
    end


    de = (r'*r)\(rnext'*rnext);
    r = rnext;
    p = rnext + p*de;

    %%% Rank check via SVD
    if rankrevtol ~= 0
            [~,D,~] = svd(p,"econ");
            d = abs(diag(D));
            rank = sum(d > rankrevtol*d(1));
        if rank ~= curr_size
            numrank = [numrank;[rank, i]];
        end
        curr_size = rank;
    end

    if coef
        ST = sigma\theta;
        tau = ga\ST;
        d = theta'*sigma*tau;
        alpha = d + ell*beta';
        try
            sigma = chol(rnext'*rnext);
        catch
            warning('breakdown, rank-deficient r, iteration %d\n',i);
            T(end-m+1:end,end-m+1:end) = alpha;
            break
        end
        [theta,beta] = qrp(sigma*tau);
        ell = theta'*sigma*ST;
        
        if i == maxit
            T = makeT(T,alpha,[]);
        else
            T = makeT(T,alpha,beta);
        end

    end
end