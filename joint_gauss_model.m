function out_warp = joint_gauss_model(out_warp, n, no)
% JOINT_GAUSS_MODEL Gaussian gnerative model from jointFPCA
% -------------------------------------------------------------------------
% This function models the functional data using a Gaussian model extracted
% from the principal components of the srvfs and warping functinos
%
% Usage: samples = joint_gauss_model(out_warp, n, no)
%
% Inputs:
% fn: matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
% time: vector of size \eqn{N} describing the sample points
% qn: matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned srvfs
% gam: warping functions
% n: number of random samples (n = 1)
%
% Output:
% Structure containing
% fs: random aligned samples
% gams: random warping function samples
% ft: random function samples

%% Separated and Warped Data
% sampling from the estimated model
fn = out_warp.fn;
time = out_warp.time;
qn = out_warp.qn;
[M, ~] = size(fn);

%% Perform PCA
jfpca = fdajpca(out_warp);
jfpca = jfpca.calc_fpca(no);
s = jfpca.latent;
U = jfpca.U;
C = jfpca.C;
mu_psi = jfpca.mu_psi;

mq_new = mean(qn,2);
id = jfpca.id;
m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));  % scaled version
mqn = [mq_new; mean(m_new)];

%% Generate Random samples
vals = mvnrnd(zeros(1,length(s)),diag(s), n);

tmp = U * vals.';
qhat = repmat(mqn, 1, n) + tmp(1:(M+1),:);
tmp = U * (vals.'/C);
vechat = tmp((M+2):end,:);
psihat = zeros(M, n);
gamhat = zeros(M, n);
for ii = 1:n
    psihat(:,ii) = exp_map(mu_psi,vechat(:,ii));
    gam_tmp = cumtrapz(linspace(0,1,M), psihat(:,ii).*psihat(:,ii));
    gamhat(:,ii) = (gam_tmp - min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
end

ft = zeros(M,n);
fhat = zeros(M,n);
for ii = 1:n
    fhat(:,ii) = cumtrapzmid(time, qhat(1:M,ii).*abs(qhat(1:M,ii)), sign(qhat(M+1,ii)).*(qhat(M+1,ii)^2), id);
    ft(:,ii) = warp_f_gamma(fhat(:,ii), gamhat(:,ii), linspace(0,1,M));
end

out_warp.fs = fhat;
out_warp.gams = gamhat;
out_warp.ft = ft;
out_warp.qs = qhat(1:M,:);
out_warp.rsamps = true;
