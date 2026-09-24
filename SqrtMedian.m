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
d = zeros(1,N);
for k = 1:N
    psi(:,k) = sqrt(max(gradient(gam(:,k),binsize),0));
    [v(:,k), d(k)] = inv_exp_map(psi_median,psi(:,k));
end
vbar = weiszfeld_step(v,d);
vbar_norm(r) = L2norm(vbar);

% compute phase median by iterative algorithm
while (vbar_norm(r) > 0.00000001 && r<501)
    psi_median = exp_map(psi_median, stp*vbar);
    r = r + 1;
    for k = 1:N
        [v(:,k), d(k)] = inv_exp_map(psi_median,psi(:,k));
    end
    vbar = weiszfeld_step(v,d);
    vbar_norm(r) = L2norm(vbar);
end

vec = v;
gam_median = cumtrapz(t,psi_median.^2)';
gam_median = (gam_median-min(gam_median))/(max(gam_median)-min(gam_median));

end

function vbar = weiszfeld_step(v,d)
% inverse-distance weighted mean of shooting vectors; observations that
% coincide with the current median carry no direction and are skipped
idx = d > 1e-10;
if ~any(idx)
    vbar = zeros(size(v,1),1);
    return
end
vbar = sum(v(:,idx)./d(idx),2)/sum(1./d(idx));
end
