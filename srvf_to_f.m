function f = srvf_to_f(q,time,fo)
% Convert SRSF to f
%
% This function converts SRSFs to functions
%
% @param q matrix of srsf
% @param time time
% @param f0 initial value of f
% @return f matrix of functions
integrand = q.*abs(q);
f = cumtrapz(time, integrand);
for i = 1:length(fo)
    f(:,i) = fo(i) + f(:,i);
end
