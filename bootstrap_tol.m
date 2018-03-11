function bounds = bootstrap_tol(out_warp, quants, nboot, alpha)
% This function models the functional data using a Gaussian model extracted from
% the principal components of the srvfs
%
% Inputs
% time vector of size \eqn{N} describing the sample points
% fn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
% qn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned srvfs
% gam warping functions
% quant array of quantiles of normal - example [.0275, 0.975]
% nboot number of bootstraps
% alpha significance level - example .05
% sort_samples sort samples (default = F)
%
% Return
% Structure containing
% \item{lwrtol}{lower tolerance fs}
% \item{uprtol}{upper tolerance fs}
% \item{mn}{tolerance of mean fs}
% \item{lwrtol_gam}{lower tolerance gams}
% \item{uprtol_gam}{upper tolerance gams}
% \item{mn_gam}{tolerance of mean gams}

time = out_warp.time;
fn = out_warp.fn;
qn = out_warp.qn;
gam = out_warp.gam;

% bootstrap CI's
bootlwr = zeros(length(time),nboot);
bootupr = zeros(length(time),nboot);
bootmean = zeros(length(time),nboot);

bootlwr_gam = zeros(length(time),nboot);
bootupr_gam = zeros(length(time),nboot);
bootmean_gam = zeros(length(time),nboot);
for ii = 1:nboot
    if (mod(ii,100)==0)
        display(ii)
    end
    samples = gauss_model(fn, time, qn, gam, size(fn,2), false);
    
    bootlwr(:,ii) = mean(samples.fs,2)+norminv(quants(1),0,1)*sqrt(var(samples.fs,0,2));
    bootupr(:,ii) = mean(samples.fs,2)+norminv(quants(2),0,1)*sqrt(var(samples.fs,0,2));
    bootmean(:,ii) = mean(samples.fs,2);
    
    bootlwr_gam(:,ii) = mean(samples.gams)+norminv(quants(1),0,1)*sqrt(var(samples.gams));
    bootupr_gam(:,ii) = mean(samples.gams)+norminv(quants(2),0,1)*sqrt(var(samples.gams));
    bootmean_gam(:,ii) = mean(samples.gams);
    
end

% get CI's from bootstraps
bounds.lwrtol = quantile(bootlwr,alpha/2,2);
bounds.uprtol = quantile(bootupr,1-alpha/2,2);
bounds.mnCI = quantile(bootmean,[alpha/2, 1-alpha/2],2);

bounds.lwrtol_gam = quantile(bootlwr_gam,alpha/2,2);
bounds.uprtol_gam = quantile(bootupr_gam,1-alpha/2,2);
bounds.mnCI_gam = quantile(bootmean_gam,[alpha/2, 1-alpha/2],2);
