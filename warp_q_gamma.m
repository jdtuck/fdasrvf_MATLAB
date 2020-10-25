function qo = warp_q_gamma(q,gamma,t,spline)
% WARP_Q_GAMMA Warp SRVF by gamma
% -------------------------------------------------------------------------
% This function warps SRVF \eqn{q} by \eqn{\gamma}
%
% Usage: qo = warp_q_gamma(q,gamma,t,spline)
%   
% Input:
% q: vector SRVF
% gamma: vector warping function
% t: time
% spline use spline interpolation (default false)
%
% Output
% qo: warped SRVF
M = length(gamma);
gam_dev = gradient(gamma, 1/(M-1));

if nargin < 4
    spline = false;
end

if spline
    qo = interp1(t, q, (t(end)-t(1)).*gamma + t(1), 'makima')'.*sqrt(gam_dev');
else
    qo = interp1(t, q, (t(end)-t(1)).*gamma + t(1))'.*sqrt(gam_dev');
end
