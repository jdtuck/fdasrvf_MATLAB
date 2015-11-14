function qo = warp_q_gamma(q,gamma,t)
M = length(gamma);
gam_dev = gradient(gamI, 1/(M-1));
qo = interp1(t, q, (t(end)-t(1)).*gamma + t(1))'.*sqrt(gam_dev');