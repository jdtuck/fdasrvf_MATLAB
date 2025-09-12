function [resmat, PNS] = PNS_warping(gam)

[d,n] = size(gam);
t = linspace(0, 1, d);
psi = zeros(d,n);
binsize = mean(diff(t));
for i=1:n
    psi(:,i) = sqrt(gradient(gam(:,i),binsize));
end
radius = mean(sqrt(sum(psi.^2)));
pnsdat = psi./repmat(sqrt(sum(psi.^2)),d,1);

[resmat, PNS]=PNSmainHDLSS(pnsdat,1);
PNS.radius = radius;