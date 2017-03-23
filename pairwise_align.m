function [f2_align, gam] = pairwise_align(f1,f2,time)
% aligns f2 to f1
q1 = f_to_srvf(f1,time);
q2 = f_to_srvf(f2,time);

gam = optimum_reparam(q1,q2,time,0,'DP',0.0,0.0,0.0);

f2_align = warp_f_gamma(f2,gam,time);
