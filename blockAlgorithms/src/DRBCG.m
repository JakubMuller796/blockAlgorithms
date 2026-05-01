function [x,T,dr_omega,numrank,conv] = DRBCG(A,b,x_0,coef,reo,tol,rankrevtol,maxit,type,xex)

tic
[n,m] = size(b);
dr_omega = [];
numrank = [];
den = norm(b,'fro');
T = [];
conv = [0,0];

x = x_0;
r = b - A*x_0;

if m == 1
    % Vector
    sigm = norm(r);
    w = r/sigm;
elseif type == "qr"
    % QR
    [w,sigm] = qrp(r);
elseif type == "svd"
    % SVD
    [aa,bb,cc] = svd(r,"econ");
    w = aa;
    sigm = bb * cc';
end

numrank = [m,0];
curr_size = m;
s = w;

if reo, V = w; end  

if coef
    beta = sigm;
    theta = eye(m); 
    ell = zeros(m,m);
end


for i = 1:maxit
    ksi = inv(s'*A*s);
    
    x = x + s*ksi*sigm;
    
    dr_omega = [dr_omega trace((x-xex)'*A*(x-xex))];
    
    if tol > 0 && conv(1) == 0
        r = w*sigm;
        res = norm(r,'fro') / den;
        if res < tol
            % fprintf("convergence in %d\n",i)
            conv(1) = toc;
            conv(2) = i;
        end
    end

    wAs = w - (A*s)*ksi;

    if reo, [wAs,V] = reorth(wAs,V,reo); end 

    if m == 1
        % Vector
        zeta = norm(wAs);
        w = (wAs)/zeta;

    elseif type == "qr"
        % QR
        [w,zeta] = qrp(wAs);

        if rankrevtol > 0
            [~,D,~] = svd(sigm,"econ");
            d = abs(diag(D));
            rank = sum(d > rankrevtol*d(1));
            if rank ~= curr_size
                numrank = [numrank;[rank, i]];
                curr_size = rank;
                % fprintf("rank-deficiency in %d\n",i)
            end
        end
    elseif type == "svd"
        % SVD
        [aa,bb,cc] = svd(wAs,"econ");
        w = aa;
        zeta = bb * cc';

        if rankrevtol > 0
            [~,D,~] = svd(sigm,"econ");
            d = abs(diag(D));
            rank = sum(d > rankrevtol*d(1));
            if rank ~= curr_size
                numrank = [numrank;[rank, i]];
                curr_size = rank;
                % fprintf("rank-deficiency in %d\n",i)
            end
        end
    end

    
    if coef
        tt = (s'*(A*s))*theta;
        d = theta'*tt;
        alpha = d + ell*beta';
        ell = zeta*theta;
        [theta,beta] = qrp(zeta*tt);
        ell = theta'*ell;
        if i == maxit
            T = makeT(T,alpha,[]);
        else
            T = makeT(T,alpha,beta);
        end
    end

    s = w + s*zeta';
    sigm = zeta*sigm;

end

end