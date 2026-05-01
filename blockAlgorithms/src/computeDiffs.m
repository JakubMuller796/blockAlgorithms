function [diffs] =  computeDiffs(T_lan, T_cg, m)

diffs = [];
n = min(size(T_lan,1),size(T_cg,1));
K = n / m;

for k = 1:K
    idx = 1:k*m;
    
    T_lan_curr = T_lan(idx,idx);
    T_cg_curr = T_cg(idx,idx);
    
    if k > 1 && k < K
        cut = (k-1)*m + (1:m);
        T_lan_curr(cut,cut) = 0;
        T_cg_curr(cut,cut) = 0;
        den = norm(T_lan_curr,"fro");
        D = T_lan_curr - T_cg_curr;
    else
        D = T_lan_curr - T_cg_curr;
        den = norm(T_lan_curr,"fro");
    end
    
    val = norm(D, 'fro');
    
    diffs = [diffs,val/den];
end

end