function f = srvf_to_f(q,time,fo)
% SRVF_TO_F Convert SRSF to f
% -------------------------------------------------------------------------
% This function converts SRSFs to functions
%
% Usage: f = srvf_to_f(q,time,fo)
%
% Input:
% q: matrix of srsf
% time: time
% fo: initial value of f (one per function)
%
% Output:
% f: matrix of functions
[M, N] = size(q);
f = zeros(M,N);
tt = time(:);
for i = 1:N
    integrand = q(:,i).*abs(q(:,i));
    % integrate the interpolating cubic spline of the derivative; exact
    % inverse of the spline derivative used in f_to_srvf (smooth=false)
    F = fnval(fnint(csapi(tt, integrand)), tt);
    f(:,i) = fo(i) + F - F(1);
end
