function fo = warp_f_gamma(f,gamma,t,spline)
if nargin < 4
    spline = false;
end

if spline
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1), 'maxima')';
else
    fo = interp1(t, f, (t(end)-t(1)).*gamma + t(1))';
end