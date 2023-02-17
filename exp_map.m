function [expgam] = exp_map(psi, v)
% EXP_MAP Exponential Map
psi = psi(:);
v = v(:);
v_norm = L2norm(v);
if sum(v_norm) == 0
    expgam = cos(v_norm) * psi;
else
    expgam = cos(v_norm) * psi + sin(v_norm) * v / v_norm;
end
