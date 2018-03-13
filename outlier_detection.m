function q_outlier = outlier_detection(q, time, mq, k)
% OUTLIER_DETECTION Outlier Detection
% -------------------------------------------------------------------------
% This function calculates outlier's using geodesic distances of the SRVFs from
% the median
%
% Usage: [mu,gam_mu,psi,vec] = outlier_detection(q, time, mq)
%        [mu,gam_mu,psi,vec] = outlier_detection(q, time, mq, k)
%
% Input:
% q: matrix (\eqn{N} x \eqn{M}) of \eqn{M} SRVF functions with \eqn{N} samples
% time: vector of size \eqn{N} describing the sample points
% mq: median calculated using time_warping
% k: cutoff threshold (default = 1.5)
%
% Output
% q_outlier: outlier functions
if nargin < 4
  k = 1.5;
end
N = size(q,2);
ds = zeros(1,N);
for kk = 1:N
  ds(kk) = sqrt(sum(trapz(time, (mq-q(:,kk)).^2)));
end

quartile_range = quantile(ds, [0 .25 .50 .75 1]);
IQR = quartile_range(4) - quartile_range(2);

thresh = quartile_range(4) + k * IQR;

ind = ds > thresh;

q_outlier = q(:,ind);
