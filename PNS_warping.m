function [resmat, PNS] = PNS_warping(gam, fast)

arguments
    gam double
    fast=true
end

[d,n] = size(gam);
t = linspace(0, 1, d);
psi = zeros(d,n);
binsize = mean(diff(t));
for i=1:n
    psi(:,i) = sqrt(gradient(gam(:,i),binsize));
end
radius = mean(sqrt(sum(psi.^2)));
pnsdat = psi./repmat(sqrt(sum(psi.^2)),d,1);

if fast
    [resmat, PNS] = fastpns(pnsdat, 1);
else
    [resmat, PNS] = PNSmainHDLSS(pnsdat,1);
end
PNS.radius = radius;
