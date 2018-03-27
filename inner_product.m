function [ip] = inner_product(psi1, psi2)
% INNER_PRODUCT Functional inner product
psi1 = psi1(:);
psi2 = psi2(:);
M = length(psi1);
t = linspace(0,1,M);
t = t(:);
ip = trapz(t,psi1.*psi2);
