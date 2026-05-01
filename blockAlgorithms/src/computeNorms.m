function [norms] =  computeNorms(omega,den)

norms = [];
for i = 1:size(omega,2)
    norms = [norms sqrt(omega(i)/den)];
end
norms = real(norms);