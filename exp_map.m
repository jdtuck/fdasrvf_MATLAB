function [expgam] = exp_map(psi, v)

v_norm = L2norm(v);
expgam = cos(v_norm) * psi + sin(v_norm) * v / v_norm;
