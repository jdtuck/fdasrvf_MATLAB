function [ip] = inner_product(psi1, psi2)

M = length(psi1);
t = linspace(0,1,M);
psi1 = psi1(:);
psi2 = psi2(:);
t = t(:);
ip = trapz(t,psi1.*psi2);
