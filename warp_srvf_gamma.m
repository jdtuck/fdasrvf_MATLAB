function qn = warp_srvf_gamma(q,gamma,scale,spline)
% WARP_SRVF_GAMMA Warp srvf by gamma
% -------------------------------------------------------------------------
% This function warps srvf \eqn{q} by \eqn{\gamma}
%
% Usage: warp_curve_gamma = warp_srvf_gamma(f,gamma,spline)
%
% Input:
% q: matrix (n,T) defining T points on n dimensional srvf
% gamma: vector warping function
% scale: scale to unit length
% spline use spline interpolation (default false)
%
% Output
% qn: warped srvf

arguments
    q
    gamma
    scale = true;
    spline = false;
end

[n,T] = size(q);
qn = q;
gam_dev = gradient(gamma, 1/T);

if spline
    for j=1:n
        qn(j,:) = interp1(linspace(0,1,T) , q(j,:), gamma, 'makima').*sqrt(gam_dev);
    end
else
    for j=1:n
        qn(j,:) = interp1(linspace(0,1,T) , q(j,:), gamma, 'linear').*sqrt(gam_dev);
    end
end

lenq = sqrt(InnerProd_Q(qn,qn));

if scale
    qn = qn/lenq;
end
