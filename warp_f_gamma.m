function fo = warp_f_gamma(f,gamma,t,spline)
% WARP_F_GAMMA Warp Function by gamma
% -------------------------------------------------------------------------
% This function warps function \eqn{f} by \eqn{\gamma}
%
% Usage: fo = warp_f_gamma(f,gamma,t,spline)
%   
% Input:
% f: vector function
% gamma: vector warping function
% t: time
% spline use spline interpolation (default false)
%
% Output
% fo: warped function
if nargin < 4
    spline = false;
end

if spline
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1), 'maxima')';
else
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1))';
end
