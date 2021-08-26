function [dy, dx] = elastic_distance(f1, f2, time, lambda, method)
% ELASTIC_DISTANCE Calculates the two elastic distances between two
% functions
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
% method: controls which optimization method (default="DP") options are
% Dynamic Programming ("DP"), Coordinate Descent ("DP2"), and Riemannian BFGS
% ("RBFGSM")
%
% Output
% dy: amplitude distance
% dx: phase distance
if nargin < 4
    lambda = 0;
    method = 'DP';
elseif nargin < 5
    method = 'DP';
end

q1 = f_to_srvf(f1,time);
q2 = f_to_srvf(f2,time);
gam = optimum_reparam(q1,q2,time.',lambda,method);
fw = warp_f_gamma(f2,gam(:),time.');
qw = f_to_srvf(fw.',time);
dy = sqrt(trapz(time,(q1-qw).^2));

time1 = linspace(0,1,length(time));
binsize = mean(diff(time1));
psi = sqrt(gradient(gam,binsize));
dx = acos(trapz(time1,psi));
