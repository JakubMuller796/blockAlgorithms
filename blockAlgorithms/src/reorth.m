function [r,V] = reorth(r,V,reo)


for k = 1:reo 
    r = r - V * (V'*r);
end
[W,~] = qr(r,"econ");
V = [V,W];
