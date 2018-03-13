function qo = warp_q_gamma(q,gamma,t,spline)
% Warp SRSF
%
% This function warps srsf \eqn{q} by \eqn{\gamma}
%
% @param q vector
% @param time time
% @param gamma vector warping function
% @param spl.int use spline interpolation (default F)
% @return qnew warped function
M = length(gamma);
gam_dev = gradient(gamma, 1/(M-1));

if nargin < 4
    spline = false;
end

if spline
    qo = interp1(t, q, (t(end)-t(1)).*gamma + t(1), 'maxima')'.*sqrt(gam_dev');
else
    qo = interp1(t, q, (t(end)-t(1)).*gamma + t(1))'.*sqrt(gam_dev');
end
