function [f2_align, gam] = pairwise_align(f1,f2,time)
% Align two functions
%
% This function aligns two functions using SRSF framework. It will align f2
% to f1
% input:
% f1: vecotr defining M samples of function 1
% f2: vecotr defining M samples of function 2
% time : time vector of length M
%
% output:
% f2_align % aligned f2
% gam % warping function
q1 = f_to_srvf(f1,time);
q2 = f_to_srvf(f2,time);

gam = optimum_reparam(q1,q2,time,0,'DP1',0.0,0.0,0.0);

f2_align = warp_f_gamma(f2,gam,time);
