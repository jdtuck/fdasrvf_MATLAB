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
% f0: initial value of f
%
% Output:
% f: matrix of functions
[M, N] = size(q);
f = zeros(M,N);
for i = 1:N
    integrand = q(:,i).*abs(q(:,i));
    f(:,i) = cumtrapz(time, integrand);
    f(:,i) = fo(i) + f(:,i);
end

