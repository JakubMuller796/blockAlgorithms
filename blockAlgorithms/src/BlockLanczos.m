function [V,T,numrank] = BlockLanczos(A, v, maxit, coef, reo, rankrevtol)

[n,m] = size(v);
T = [];

V_prev = zeros(n,m);
[V,beta] = qrp(v);
V_k = V;

numrank = [size(V_k,2), 0];
curr_size = size(V_k,2);

for k = 1:maxit

    w = A*V_k - V_prev*beta';
    alpha = V_k'*w;
    w = w - V_k*alpha;
    
    if reo > 0
        w_reo = w;
        for r = 1:reo
            w_reo = w_reo - V*(V'*w_reo);
        end
        w = w_reo;
    end
        
    % rank monitor
    if rankrevtol > 0
            [~,D,~] = svd(w,"econ");
            d = abs(diag(D));
            rank = sum(d > rankrevtol*d(1));
        if rank ~= curr_size
            numrank = [numrank;[rank, k]];
        end
        curr_size = rank;
    end

    [V_next,beta] = qrp(w);
    
    % save to T
    if k == maxit && coef == 1
        T = makeT(T,alpha,[]);
    elseif coef == 1
        T = makeT(T,alpha,beta);
    end
    
    if coef
        V = [V V_next];
    end
    
    V_prev = V_k;
    V_k = V_next;
end

end