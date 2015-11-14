function f = srvf_to_f(q,time,fo)
integrand = q.*abs(q);
f = cumtrapz(time, integrand);
for i = 1:length(fo)
    f(:,i) = fo(i) + f(:,i);
end