function [gam_median, psi_median, psi, vec] = SqrtMedian(gam)
% SQRTMEDIAN SRVF transform of warping functions calculates median
% -------------------------------------------------------------------------
% This function calculates the srvf of warping functions with corresponding
% shooting vectors and finds the median
%
% Usage: [gam_median, psi_median, psi] = SqrtMedian(gam)
%
% Input
% gam: matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
% 
% Output:
% median: Karcher median psi function
% gam_median: Karcher mean warping function
% psi: srvf of warping functions
% vec: shooting vectors
[M, N] = size(gam);
t = linspace(0,1,M);

% Initialization
psi_median(1:M,1) = 1;
r = 1; stp = 0.3;

% compute psi-functions
binsize = mean(diff(t));
psi = zeros(M,N);
v = zeros(M,N);
vtil = v;
d = zeros(1,N);
dtil = zeros(1,N);
for k = 1:N
    psi(:,k) = sqrt(gradient(gam(:,k),binsize));
    [v(:,k), d(k)] = inv_exp_map(psi_median,psi(:,k));
    vtil(:,k)=v(:,k)/d(k);
    dtil(k)=1/d(k);
end
vbar=sum(vtil,2)*sum(dtil)^(-1);
vbar_norm(r) = L2norm(vbar);

% compute phase median by iterative algorithm
while (vbar_norm(r) > 0.00000001 && r<501)
    psi_median = exp_map(psi_median, stp*vbar);
    r = r + 1;
    for k = 1:N
        v(:,k) = inv_exp_map(psi_median,psi(:,k));
        d(k) = acos(inner_product(psi_median,psi(:,k)));
        vtil(:,k)=v(:,k)/d(k);
        dtil(k)=1/d(k);
    end
    vbar=sum(vtil,2)*sum(dtil)^(-1);
    vbar_norm(r) = L2norm(vbar);
end

vec = v;
gam_median = cumtrapz(t,psi_median.^2)';
gam_median = (gam_median-min(gam_median))/(max(gam_median)-min(gam_median));
