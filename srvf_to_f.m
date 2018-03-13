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
integrand = q.*abs(q);
f = cumtrapz(time, integrand);
for i = 1:length(fo)
    f(:,i) = fo(i) + f(:,i);
end
