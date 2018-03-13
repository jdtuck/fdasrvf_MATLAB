function [dy, dx] = elastic_distance(f1, f2, time, lambda)
% ELASTIC_DISTANCE Calculates two elastic distances
% -------------------------------------------------------------------------
% This functions calculates the distances between functions,
% \eqn{D_y} and \eqn{D_x}, where function 1 is aligned to function 2
%
% Usage: [dy, dx] = elastic_distance(f1, f2, time)
%        [dy, dx] = elastic_distance(f1, f2, time, lambda)
%
% Input:
% f1: sample function 1
% f2: sample function 1
% time: vector of size \eqn{N} describing the sample points
% lambda controls amount of warping (default = 0)
%
% Output
% dy: amplitude distance
% dx: phase distance
if nargin < 4
  lambda = 0;
end
q1 = f_to_srvf(f1,time);
q2 = f_to_srvf(f2,time);
gam = optimum_reparam(q1,q2,time,lambda);
fw = warp_f_gamma(f2,gam,time);
qw = f_to_srvf(fw,time);
dy = sqrt(trapz(time,(q1-qw).^2));

time1 = linspace(0,1,length(time));
binsize = mean(diff(time1));
psi = sqrt(gradient(gam,binsize));
v = inv_exp_map(ones(1,length(gam)), psi);
dx = sqrt(trapz(time,v.^2));
