function fo = warp_f_gamma(f,gamma,t,spline)
% Warp Function
%
% This function warps function \eqn{f} by \eqn{\gamma}
%
% @param f vector function
% @param time time
% @param gamma vector warping function
% @param spl.int use spline interpolation (default F)
% @return fnew warped function
if nargin < 4
    spline = false;
end

if spline
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1), 'maxima')';
else
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1))';
end
