function bounds = bootstrap_tol(out_warp, quants, nboot, alpha)
% BOOTSTRAP_TOL Boostrapped tolerance bounds
% -------------------------------------------------------------------------
% 
% This function models the functional data using a Gaussian model extracted 
% from the principal components of the srvfs
% 
% Usage:  bounds = bootstrap_tol(out_warp, quants, nboot, alpha)
%
% Inputs:
% out_warp: fdawarp object of aligned data
% quant: array of quantiles of normal - example [.0275, 0.975]
% nboot: number of bootstraps
% alpha: significance level - example .05
%
% Output:
% Structure containing
% lwrtol: lower tolerance fs
% uprtol: upper tolerance fs
% mnCI: tolerance of mean fs
% lwrtol_gam: lower tolerance gams
% uprtol_gam: upper tolerance gams
% mnCI_gam: tolerance of mean gams

if (~isa(out_warp,'fdawarp'))
    error('Require input of class fdawarp');
end

if (isempty(out_warp.fn))
    error('Please align using method time_warping');
end
fn = out_warp.fn;
time = out_warp.time;

% bootstrap CI's
bootlwr = zeros(length(time),nboot);
bootupr = zeros(length(time),nboot);
bootmean = zeros(length(time),nboot);

bootlwr_gam = zeros(length(time),nboot);
bootupr_gam = zeros(length(time),nboot);
bootmean_gam = zeros(length(time),nboot);
parpool();
fprintf('Bootstrap Sampling. \n');
obj = ProgressBar(nboot, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd, ...
    'Title', 'Progress' ...
    );

obj.setup([], [], []);
parfor ii = 1:nboot
    samples = gauss_model(out_warp, size(fn,2), false);
    
    bootlwr(:,ii) = mean(samples.fs,2)+norminv(quants(1),0,1)*sqrt(var(samples.fs,0,2));
    bootupr(:,ii) = mean(samples.fs,2)+norminv(quants(2),0,1)*sqrt(var(samples.fs,0,2));
    bootmean(:,ii) = mean(samples.fs,2);
    
    bootlwr_gam(:,ii) = mean(samples.gams)+norminv(quants(1),0,1)*sqrt(var(samples.gams));
    bootupr_gam(:,ii) = mean(samples.gams)+norminv(quants(2),0,1)*sqrt(var(samples.gams));
    bootmean_gam(:,ii) = mean(samples.gams);
    updateParallel([], pwd);
    
end
obj.release();

delete(gcp('nocreate'))

% get CI's from bootstraps
bounds.lwrtol = quantile(bootlwr,alpha/2,2);
bounds.uprtol = quantile(bootupr,1-alpha/2,2);
bounds.mnCI = quantile(bootmean,[alpha/2, 1-alpha/2],2);

bounds.lwrtol_gam = quantile(bootlwr_gam,alpha/2,2);
bounds.uprtol_gam = quantile(bootupr_gam,1-alpha/2,2);
bounds.mnCI_gam = quantile(bootmean_gam,[alpha/2, 1-alpha/2],2);
