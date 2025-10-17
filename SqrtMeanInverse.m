function gamI = SqrtMeanInverse(gam)
% SQRTMEANINVERSE SRVF transform of warping functions calculates mean
% inverse
% -------------------------------------------------------------------------
% This function calculates the srvf of warping functions with corresponding
% shooting vectors and finds the inverse of the mean
%
% Usage: gamI = SqrtMeanInverse(gam)
%
% Input:
% gam: matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
%
% Output
% gamI: Inverse Karcher mean warping function

[T,n] = size(gam);
time = linspace(0,1,T);

psi = zeros(T,n);
binsize = mean(diff(time));
for i=1:n
    psi(:,i) = sqrt(gradient(gam(:,i),binsize));
end

%Find direction
mnpsi = mean(psi,2);
dqq = sqrt(sum((psi - mnpsi*ones(1,n)).^2,1));
[~, min_ind] = min(dqq);
mu = psi(:,min_ind);
maxiter = 501;
lvm = zeros(1,maxiter);
vec = zeros(T,n);
stp = .3;
iter = 1;

for i = 1:n
    vec(:,i) = inv_exp_map(mu,psi(:,i));
end
vbar=mean(vec,2);
lvm(iter) = L2norm(vbar);

while (lvm(iter) > 0.00000001 && iter<maxiter)
    mu = exp_map(mu, stp*vbar);
    iter = iter + 1;
    for i = 1:n
        vec(:,i) = inv_exp_map(mu,psi(:,i));
    end
    vbar=mean(vec,2);
    lvm(iter) = L2norm(vbar);
end

gam_mu = cumtrapz(time,mu.^2);

gam_mu = (gam_mu-min(gam_mu))/(max(gam_mu)-min(gam_mu));
gamI = invertGamma(gam_mu);
