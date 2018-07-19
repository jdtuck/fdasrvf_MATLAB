function fn = warp_curve_gamma(f,gamma,spline)
% warp_curve_gamma Warp curve by gamma
% -------------------------------------------------------------------------
% This function warps curve \eqn{f} by \eqn{\gamma}
%
% Usage: warp_curve_gamma = warp_curve_gamma(f,gamma,spline)
%
% Input:
% f: matrix (n,T) defining T points on n dimensional curve
% gamma: vector warping function
% spline use spline interpolation (default false)
%
% Output
% fn: warped function
if nargin < 4
    spline = false;
end

fn = f;
if spline
    for j=1:n
        fn(j,:) = interp1(linspace(0,1,T) , f(j,:) ,gamma,'spline');
    end
else
    for j=1:n
        fn(j,:) = interp1(linspace(0,1,T) , f(j,:) ,gamma,'linear');
    end
end
