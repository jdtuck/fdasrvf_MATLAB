function fo = warp_f_gamma(f,gamma,t)
fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1))';