function [out, chains] = pairwise_align_bayes_infHMCbasis(y1i, y2i, time, mcmcopts)
% PAIRWISE_ALIGN_BAYES_INFHMCBASIS Align two trajectories using hierarchical
% Bayesian framework assuming measuring error
% -------------------------------------------------------------------------
%
% Usage:  out = pairwise_align_bayes_infHMCbasis(y1i, y2i, time)
%         out = pairwise_align_bayes_infHMCbasis(y1i, y2i, time, mcmcopts)
%
% Input:
% y1i: vector defining M samples of function 1
% y2i: vector defining M samples of function 2
% time: time vector of length M
%
% default mcmc options
% mcmcopts.iter = 1e4;
% mcmcopts.nchains = 1;
% mcmcopts.burnin = min(5e3,mcmcopts.iter/2); %must be >1
% mcmcopts.alpha0 = 0.1;
% mcmcopts.beta0 = 0.1;
% mcmcopts.alpha = 1;
% mcmcopts.beta = 1;
% mcmcopts.h = .01; % HMC step size is sqrt(h) (rotation angle)
% mcmcopts.L = 4;  % maximum number of leapfrog steps
% mcmcopts.adapt = true; % adapt the step size during burn-in
% mcmcopts.target_accept = 0.65; % target HMC acceptance rate when adapting
% mcmcopts.f1propvar = 0.0001;
% mcmcopts.f2propvar = 0.0001; % proposal variance for f2
% mcmcopts.L1propvar = 0.3; %controls L1 (in f1 cov) acceptance rate
% mcmcopts.L2propvar = 0.3;  %controls L2 (in f2 cov) acceptance rate
% mcmcopts.vpriorvar = 1; % prior variance scale of the warping coefficients
% mcmcopts.vpriordecay = 1; % prior variance of pair j is vpriorvar/j^vpriordecay
% mcmcopts.initcoef = repelem(0, 20).'; % init coef
% mcmcopts.basis = 'fourier';  %'fourier' 'legendre'
% mcmcopts.warp_rep = 'shooting'; % 'shooting' (v-transform) or 'clr' (h-transform)
% mcmcopts.extrainfo = true;
% mcmcopts.parallel = true; % run chains in parallel when nchains > 1
% mcmcopts.nworkers = 0; % number of workers for parallel pool (0: min(numcores, nchains))
%
% Output:
% structure containing
% out.f2_warped: aligned f2
% out.gamma: warping function
% out.v_coef: final v_coef
% out.psi: final psi
% out.sigma1: final sigma
%
% if extrainfo
% out.theta_accept: accept of psi samples
% out.f2_accept: accept of f2 samples
% out.SSE: SSE of level 2 likelihood
% out.gamma_mat: posterior gammas
% out.gamma_stats: posterior gamma stats
% out.phasedist: phase distance posterior
% out.ampdist: amplitude distance posterior
% out.f2_warped: f2 aligned functions
% out.epsilon: HMC step size (rotation angle) at each iteration
%
% Warping functions: with warp_rep = 'shooting', v = Phi*c is a shooting
% vector at the identity of the Hilbert sphere of SRSFs of warps and must
% satisfy |v| < pi (draws outside are rejected); with warp_rep = 'clr',
% h = Phi*c is the centered log-derivative, gam' = exp(h)/int exp(h), which
% is unconstrained. The likelihood gradient with respect to c is exact (see
% bayes_warp_geom) and v is updated with the preconditioned HMC of Beskos et
% al. (2011) (see infhmc_step). time is assumed to be uniformly sampled.
%
% J. D. Tucker, L. Shand, and K. Chowdhary. “Multimodal Bayesian Registration of Noisy Functions using 
%       Hamiltonian Monte Carlo”, Computational Statistics and Data Analysis, accepted, 2021.
%

arguments
    y1i
    y2i
    time
    mcmcopts.iter = 1e4;
    mcmcopts.burnin = []; % default min(5e3, iter/2)
    mcmcopts.thin = 50;
    mcmcopts.nchains = 1;
    mcmcopts.alpha0 = 0.1;
    mcmcopts.beta0 = 0.1;
    mcmcopts.alpha = 1;
    mcmcopts.beta = 1;
    mcmcopts.h = .01; % step size is sqrt(h)
    mcmcopts.L = 4;  % maximum number of leapfrog steps
    mcmcopts.adapt = true;
    mcmcopts.target_accept = 0.65;
    mcmcopts.f1propvar = 0.0001;
    mcmcopts.f2propvar = 0.0001; % proposal variance for f2
    mcmcopts.L1propvar = 0.3; %controls L1 (in f1 cov) acceptance rate
    mcmcopts.L2propvar = 0.3;  %controls L2 (in f2 cov) acceptance rate
    mcmcopts.vpriorvar = 1; % prior variance scale of theta
    mcmcopts.vpriordecay = 1; % prior variance decays as 1/j^vpriordecay
    mcmcopts.initcoef = repelem(0, 20).'; % init coef
    mcmcopts.nbasis = 10; %#ok<*INUSA> unused, kept for backward compatibility
    mcmcopts.basis = 'fourier';  %'fourier' 'legendre'
    mcmcopts.warp_rep {mustBeMember(mcmcopts.warp_rep,{'shooting','clr'})} = 'shooting';
    mcmcopts.sampfreq = 5; %controls smoothness of fitrgp
    mcmcopts.seed = randi(floor(intmax/10));
    mcmcopts.extrainfo = true;
    mcmcopts.parallel = true;
    mcmcopts.nworkers = 0;
end

if (length(y1i) ~= length(y2i))
    error('Length of y1 and y2 must be equal')
end
if (length(y1i) ~= length(time))
    error('Length of y1 and time must be equal')
end
if (mod(length(mcmcopts.initcoef), 2) ~= 0)
    error('Length of mcmcopts.initcoef must be even')
end
if isempty(mcmcopts.burnin)
    mcmcopts.burnin = min(5e3, floor(mcmcopts.iter/2));
end
if mcmcopts.burnin < 1 || mcmcopts.burnin >= mcmcopts.iter
    error('mcmcopts.burnin must be between 1 and mcmcopts.iter-1')
end

% setup random start points for more than 1 chain
random_starts = zeros(length(mcmcopts.initcoef),mcmcopts.nchains);
if mcmcopts.nchains > 1
    rng(mcmcopts.seed)
    for i = 1:mcmcopts.nchains
        randcoef = -1 + (2)*rand(length(mcmcopts.initcoef),1); %.1*(2*rand(length(mcmcopts.initcoef),1) - 1.0);
        if strcmpi(mcmcopts.warp_rep,'shooting') && norm(randcoef) >= 0.9*pi
            % keep the shooting vector inside the injectivity radius
            randcoef = randcoef * (0.9*pi/norm(randcoef)) * rand;
        end
        random_starts(:,i) = randcoef;
    end
end

isparallel = mcmcopts.parallel;
try
    tmp = ver('parallel');
catch
    isparallel = false;
end

if mcmcopts.nchains == 1
    isparallel = false;
end

if isparallel
    nThreads = mcmcopts.nworkers;
    if isempty(nThreads) || nThreads < 1
        nThreads = min(feature('numcores'),mcmcopts.nchains);
    end
    start_parpool(nThreads);
    mcmcopts_p = cell(1,mcmcopts.nchains);
    for i = 1:mcmcopts.nchains
        mcmcopts.initcoef = random_starts(:,i);
        mcmcopts_p{i} = mcmcopts;
    end
end

if mcmcopts.nchains > 1 && isparallel
    % update the chain for iter-1 times
    obj = ProgressBar(mcmcopts.nchains, ...
        'IsParallel', true, ...
        'WorkerDirectory', pwd, ...
        'Title', 'Running MCMC Chains' ...
        );
else
    % update the chain for iter-1 times
    obj = ProgressBar(mcmcopts.nchains, ...
        'Title', 'Running MCMC Chains' ...
        );
end

% run chains
chains = cell(1,mcmcopts.nchains);
if isparallel
    obj.setup([], [], []);
    parfor i = 1:mcmcopts.nchains
        rng(mcmcopts.seed+i)
        chains{i} = run_mcmc(y1i, y2i, time, mcmcopts_p{i});
        updateParallel([], pwd);
    end
    obj.release();
else
    for i = 1:mcmcopts.nchains
        rng(mcmcopts.seed+i)
        mcmcopts.initcoef = random_starts(:,i);
        chains{i} = run_mcmc(y1i, y2i, time, mcmcopts);
        obj([], [], []);
    end
    obj.release();
end

% combine outputs
Nsamples = size(chains{1}.f1,1);
M = size(chains{1}.f1,2);
out.f1 = zeros(Nsamples*mcmcopts.nchains,M);
out.f2 = zeros(Nsamples*mcmcopts.nchains,M);
out.gamma = zeros(M,mcmcopts.nchains);
out.v_coef = zeros(Nsamples*mcmcopts.nchains,size(chains{1}.v_coef,2));
out.psi = zeros(M,mcmcopts.nchains);
out.sigma = zeros(1,Nsamples*mcmcopts.nchains);
out.sigma1 = zeros(1,Nsamples*mcmcopts.nchains);
out.sigma2 = zeros(1,Nsamples*mcmcopts.nchains);
out.s1 = zeros(1,Nsamples*mcmcopts.nchains);
out.s2 = zeros(1,Nsamples*mcmcopts.nchains);
out.L1 = zeros(1,Nsamples*mcmcopts.nchains);
out.L2 = zeros(1,Nsamples*mcmcopts.nchains);
out.f2_warped_mu = zeros(M,mcmcopts.nchains);
out.epsilon = zeros(mcmcopts.nchains, length(chains{1}.epsilon));

if (mcmcopts.extrainfo)
    Nsamplesa = length(chains{1}.theta_accept);
    out.theta_accept = zeros(1,Nsamplesa*mcmcopts.nchains);
    out.f1_accept = zeros(1,Nsamplesa*mcmcopts.nchains);
    out.f2_accept = zeros(1,Nsamplesa*mcmcopts.nchains);
    out.L1_accept = zeros(1,Nsamplesa*mcmcopts.nchains);
    out.L2_accept = zeros(1,Nsamplesa*mcmcopts.nchains);
    out.gamma_mat = zeros(M,Nsamples*mcmcopts.nchains);
    % @todo gamma_stats...mode finding
    out.SSE = zeros(1,(Nsamplesa+1)*mcmcopts.nchains);
    out.logl = zeros(1,(Nsamplesa+1)*mcmcopts.nchains);
    out.f2_warped = zeros(Nsamples*mcmcopts.nchains,M);
    out.phasedist = zeros(1,Nsamples*mcmcopts.nchains);
    out.ampdist = zeros(1,Nsamples*mcmcopts.nchains);
end

for i = 1:mcmcopts.nchains
    a = (i-1)*Nsamples + 1;
    b = (i)*Nsamples;
    out.f1(a:b,:) = chains{i}.f1;
    out.f2(a:b,:) = chains{i}.f2;
    out.gamma(:,i) = chains{i}.gamma;
    out.v_coef(a:b,:) = chains{i}.v_coef;
    out.psi(:,i) = chains{i}.psi;
    out.sigma(a:b) = chains{i}.sigma;
    out.sigma1(a:b) = chains{i}.sigma1;
    out.sigma2(a:b) = chains{i}.sigma2;
    out.s1(a:b) = chains{i}.s1;
    out.s2(a:b) = chains{i}.s2;
    out.L1(a:b) = chains{i}.L1;
    out.L2(a:b) = chains{i}.L2;
    out.f2_warped_mu(:,i) = chains{i}.f2_warped_mu;
    out.epsilon(i,:) = chains{i}.epsilon;

    if (mcmcopts.extrainfo)
        a1 = (i-1)*Nsamplesa + 1;
        b1 = (i)*Nsamplesa;
        out.theta_accept(a1:b1) = chains{i}.theta_accept;
        out.f1_accept(a1:b1) = chains{i}.f1_accept;
        out.f2_accept(a1:b1) = chains{i}.f2_accept;
        out.L1_accept(a1:b1) = chains{i}.L1_accept;
        out.L2_accept(a1:b1) = chains{i}.L2_accept;
        out.gamma_mat(:,a:b) = chains{i}.gamma_mat;
        a1 = (i-1)*(Nsamplesa+1) + 1;
        b1 = (i)*(Nsamplesa+1);
        out.SSE(a1:b1) = chains{i}.SSE;
        out.logl(a1:b1) = chains{i}.logl;
        out.f2_warped(a:b,:) = chains{i}.f2_warped;
        out.phasedist(a:b) = chains{i}.phasedist;
        out.ampdist(a:b) = chains{i}.ampdist;
    end
end

% finding modes
if mcmcopts.nchains > 1
    Dx = zeros(mcmcopts.nchains,mcmcopts.nchains);
    time1 = linspace(0,1,size(out.gamma,1));
    binsize = mean(diff(time1));
    for i = 1:mcmcopts.nchains
        for j = i+1:mcmcopts.nchains
            psi1 = sqrt(max(gradient(out.gamma(:,i),binsize),0));
            psi2 = sqrt(max(gradient(out.gamma(:,j),binsize),0));
            Dx(i,j) = acos(min(max(trapz(time1,psi1.*psi2),-1),1));
        end
    end
    Dx = Dx + Dx.';
    
    % cluster modes
    y = squareform(Dx);
    Z = linkage(y,'complete');
    cutoff = median(Dx(:));
    T = cluster(Z,'Cutoff',cutoff,'Criterion','distance');
    N = unique(T);
    
    % find mean and confidence region of cluster
    out.posterior_gamma_modes = zeros(M,length(N));
    out.posterior_gamma_modes_cr = zeros(M,2,length(N));
    for i = 1:length(N)
        idx = find(T == i);
        tmp = zeros(M,Nsamples*length(idx));
        for j = 1:length(idx)
            a = (j-1)*Nsamples + 1;
            b = (j)*Nsamples;
            tmp(:,a:b) = chains{idx(j)}.gamma_mat;
        end
        [~,gam_mu,~,~] = SqrtMean(tmp);
        out.posterior_gamma_modes(:,i) = gam_mu;
        out.posterior_gamma_modes_cr(:,:,i) = statsFun(tmp).';
    end
end

% thining
out.f1 = out.f1(1:mcmcopts.thin:end,:);
out.f2 = out.f2(1:mcmcopts.thin:end,:);
out.v_coef = out.v_coef(1:mcmcopts.thin:end,:);
out.sigma = out.sigma(1:mcmcopts.thin:end);
out.sigma1 = out.sigma1(1:mcmcopts.thin:end);
out.sigma2 = out.sigma2(1:mcmcopts.thin:end);
out.s1 = out.s1(1:mcmcopts.thin:end);
out.s2 = out.s2(1:mcmcopts.thin:end);
out.L1 = out.L1(1:mcmcopts.thin:end);
out.L2 = out.L2(1:mcmcopts.thin:end);

if (mcmcopts.extrainfo)
    out.theta_accept = out.theta_accept(1:mcmcopts.thin:end);
    out.f1_accept = out.f1_accept(1:mcmcopts.thin:end);
    out.f2_accept = out.f2_accept(1:mcmcopts.thin:end);
    out.L1_accept = out.L1_accept(1:mcmcopts.thin:end);
    out.L2_accept = out.L2_accept(1:mcmcopts.thin:end);
    out.gamma_mat = out.gamma_mat(:,1:mcmcopts.thin:end);
    out.SSE = out.SSE(1:mcmcopts.thin:end);
    out.logl = out.logl(1:mcmcopts.thin:end);
    out.f2_warped = out.f2_warped(1:mcmcopts.thin:end,:);
    out.phasedist = out.phasedist(1:mcmcopts.thin:end);
    out.ampdist = out.ampdist(1:mcmcopts.thin:end);
end

end

% main mcmc function
function out = run_mcmc(y1i, y2i, time, mcmcopts)

% Number of sig figs to report in gamma_mat
SIG_GAM = 13;
iter = mcmcopts.iter;

y1 = y1i(:).';
y2 = y2i(:).';
time = time(:).';
T = length(time);

% normalize time to [0,1] (uniform sampling assumed)
t01 = (time-min(time))./(max(time)-min(time));

% warping function parameter settings
valid_index = mcmcopts.burnin:iter;
ncoef = length(mcmcopts.initcoef);
rep = double(strcmpi(mcmcopts.warp_rep,'clr'));
Phi = bayes_warp_basis(t01, ncoef, mcmcopts.basis);
C = (mcmcopts.vpriorvar ./ repelem(1:ncoef/2,2).^mcmcopts.vpriordecay).';
epsilon = sqrt(mcmcopts.h);
L = mcmcopts.L;

% f1, f2 GP priors, initialized from a GP fit
Dmat = squareform(pdist(time'));
model1 = fitrgp(time(1:mcmcopts.sampfreq:end)', y1(1:mcmcopts.sampfreq:end)');
model2 = fitrgp(time(1:mcmcopts.sampfreq:end)', y2(1:mcmcopts.sampfreq:end)');

s1_curr = model1.KernelInformation.KernelParameters(2)^2;
s2_curr = model2.KernelInformation.KernelParameters(2)^2;
L1_curr = model1.KernelInformation.KernelParameters(1);
L2_curr = model2.KernelInformation.KernelParameters(1);
gp1 = gp_prior(s1_curr, L1_curr, Dmat);
gp2 = gp_prior(s2_curr, L2_curr, Dmat);

% random-walk proposals for f1, f2 (fixed covariance, factor once)
R_f1prop = chol_jitter(exp2corr(mcmcopts.f1propvar, L1_curr, Dmat));
R_f2prop = chol_jitter(exp2corr(mcmcopts.f2propvar, L2_curr, Dmat));

% result objects
result.v_coef = zeros(iter,ncoef);
result.sigma = zeros(1,iter);
result.sigma1 = zeros(1,iter);
result.sigma2 = zeros(1,iter);
result.f1 = zeros(iter,T);
result.f2 = zeros(iter,T);
result.s1 = zeros(1,iter);
result.s2 = zeros(1,iter);
result.L1 = zeros(1,iter);
result.L2 = zeros(1,iter);
result.logl = zeros(1,iter);
result.SSE = zeros(1,iter);
result.epsilon = zeros(1,iter);
result.theta_accept = zeros(1,iter);
result.f2_accept = zeros(1,iter);
result.f1_accept = zeros(1,iter);
result.L2_accept = zeros(1,iter);
result.L1_accept = zeros(1,iter);

% init
v_coef_curr = mcmcopts.initcoef(:);
sigma_curr = 1;
sigma1_curr = 0.01;
sigma2_curr = 0.01;
f1_curr = predict(model1,time').';
f2_curr = predict(model2,time').';
q1_curr = f_to_srvf(f1_curr.',t01,0);
q2_curr = f_to_srvf(f2_curr.',t01,0);

[nll, ~, SSE_curr] = bayes_warp_geom(v_coef_curr, Phi, rep, q1_curr, q2_curr, sigma_curr);
if ~isfinite(nll)
    error('Initial warping coefficients are outside the injectivity radius (|v| >= pi)')
end

result.v_coef(1,:) = v_coef_curr;
result.f1(1,:) = f1_curr;
result.f2(1,:) = f2_curr;
result.sigma(1) = sigma_curr;
result.sigma1(1) = sigma1_curr;
result.sigma2(1) = sigma2_curr;
result.s1(1) = s1_curr;
result.s2(1) = s2_curr;
result.L1(1) = L1_curr;
result.L2(1) = L2_curr;
result.SSE(1) = SSE_curr;
result.logl(1) = -nll;
result.epsilon(1) = epsilon;

if mcmcopts.nchains == 1
    obj = ProgressBar(iter-1, ...
        'Title', 'Running MCMC' ...
        );
end

for m = 2:iter

    % update f1 | v, f2
    [f1_curr, q1_curr, SSE_curr, f1_accept] = update_f(f1_curr, q1_curr, SSE_curr, y1, ...
        @(q) warp_sse(v_coef_curr, Phi, rep, q, q2_curr, sigma_curr), ...
        gp1, R_f1prop, t01, sigma_curr, sigma1_curr);

    % update f2 | v, f1
    [f2_curr, q2_curr, ~, f2_accept] = update_f(f2_curr, q2_curr, SSE_curr, y2, ...
        @(q) warp_sse(v_coef_curr, Phi, rep, q1_curr, q, sigma_curr), ...
        gp2, R_f2prop, t01, sigma_curr, sigma2_curr);

    % update v by infinite-dimensional HMC; geometry evaluated at the
    % current f1, f2 and sigma
    geom = @(c) bayes_warp_geom(c, Phi, rep, q1_curr, q2_curr, sigma_curr);
    [nll, g] = geom(v_coef_curr);
    [v_coef_curr, nll, ~, ~, theta_accept] = infhmc_step(v_coef_curr, nll, g, [], ...
        @(c) geom_aux(geom, c), C, epsilon, L);
    SSE_curr = 2 * sigma_curr * (nll - T/2*log(sigma_curr));

    % adapt step size during burn-in (Robbins-Monro on log step size)
    if mcmcopts.adapt && m <= mcmcopts.burnin
        gain = 1 / m^0.6;
        epsilon = min(exp(log(epsilon) + gain*(theta_accept - mcmcopts.target_accept)), pi/2);
    end

    % update sigma^2 (SRVF model variance)
    newshape = T/2 + mcmcopts.alpha;
    newscale = 1/2 * SSE_curr + mcmcopts.beta;
    sigma_curr = 1/gamrnd(newshape,1/newscale);

    % update sigma1^2, sigma2^2 (measurement error variances)
    newshape = T/2 + mcmcopts.alpha0;
    sigma1_curr = 1/gamrnd(newshape,1/(sum((y1 - f1_curr).^2)/2 + mcmcopts.beta0));
    sigma2_curr = 1/gamrnd(newshape,1/(sum((y2 - f2_curr).^2)/2 + mcmcopts.beta0));

    % update GP variances s1^2, s2^2
    s1_curr = 1/gamrnd(newshape,1/(quad_form(gp1.Rcorr, f1_curr)/2 + mcmcopts.beta0));
    s2_curr = 1/gamrnd(newshape,1/(quad_form(gp2.Rcorr, f2_curr)/2 + mcmcopts.beta0));
    gp1 = gp_prior(s1_curr, L1_curr, Dmat, gp1.Rcorr);
    gp2 = gp_prior(s2_curr, L2_curr, Dmat, gp2.Rcorr);

    % update GP length scales L1, L2
    [L1_curr, gp1, L1_accept] = update_length(f1_curr, gp1, s1_curr, L1_curr, mcmcopts.L1propvar, Dmat);
    [L2_curr, gp2, L2_accept] = update_length(f2_curr, gp2, s2_curr, L2_curr, mcmcopts.L2propvar, Dmat);

    % save update to results
    result.v_coef(m,:) = v_coef_curr;
    result.f1(m,:) = f1_curr;
    result.f2(m,:) = f2_curr;
    result.sigma(m) = sigma_curr;
    result.sigma1(m) = sigma1_curr;
    result.sigma2(m) = sigma2_curr;
    result.s1(m) = s1_curr;
    result.s2(m) = s2_curr;
    result.L1(m) = L1_curr;
    result.L2(m) = L2_curr;
    result.SSE(m) = SSE_curr;
    result.logl(m) = -(SSE_curr/(2*sigma_curr) + T/2*log(sigma_curr));
    result.epsilon(m) = epsilon;
    result.theta_accept(m) = theta_accept;
    result.f1_accept(m) = f1_accept;
    result.f2_accept(m) = f2_accept;
    result.L1_accept(m) = L1_accept;
    result.L2_accept(m) = L2_accept;
    if mcmcopts.nchains == 1
        obj.step([], [], []);
    end
end
if mcmcopts.nchains == 1
    obj.release();
end

% posterior warping functions
[gamma_mat, psi_mat] = bayes_coef_to_gam(result.v_coef(valid_index,:).', Phi, rep);
gamma_mat = round(gamma_mat, SIG_GAM);

% posterior mean of psi (extrinsic mean, projected back to the sphere)
psi_mean = mean(psi_mat,2);
psi_mean = psi_mean / sqrt(trapz(t01, psi_mean.^2));
gam_mean = [0; cumsum(diff(t01(:)) .* (psi_mean(1:end-1).^2 + psi_mean(2:end).^2)/2)];
gam_mean = norm_gam(gam_mean);

% return object
out.f1 = result.f1(valid_index,:);
out.f2 = result.f2(valid_index,:);
out.gamma = gam_mean;
out.v_coef = result.v_coef(valid_index,:);
out.psi = psi_mean;
out.sigma = result.sigma(valid_index);
out.sigma1 = result.sigma1(valid_index);
out.sigma2 = result.sigma2(valid_index);
out.s1 = result.s1(valid_index);
out.s2 = result.s2(valid_index);
out.L1 = result.L1(valid_index);
out.L2 = result.L2(valid_index);
out.epsilon = result.epsilon;
out.f2_warped_mu = warp_f_gamma(mean(out.f2,1), gam_mean, t01);

if (mcmcopts.extrainfo)
    nv = length(valid_index);
    Dx = real(acos(min(max(trapz(t01, psi_mat, 1), -1), 1)));
    Dy = zeros(1,nv);
    f2_warped = zeros(nv,T);
    for ii = 1:nv
        q1s = f_to_srvf(out.f1(ii,:).', t01, 0);
        q2s = f_to_srvf(out.f2(ii,:).', t01, 0);
        q2warp = warp_q_gamma(q2s, gamma_mat(:,ii), t01);
        Dy(ii) = sqrt(trapz(t01,(q1s(:)-q2warp(:)).^2));
        f2_warped(ii,:) = warp_f_gamma(out.f2(ii,:)', gamma_mat(:,ii), t01);
    end
    out.theta_accept = result.theta_accept(2:end);
    out.f1_accept = result.f1_accept(2:end);
    out.f2_accept = result.f2_accept(2:end);
    out.L1_accept = result.L1_accept(2:end);
    out.L2_accept = result.L2_accept(2:end);
    out.gamma_mat = gamma_mat;
    out.gamma_stats = statsFun(gamma_mat).';
    out.SSE = result.SSE;
    out.logl = result.logl;
    out.phasedist = Dx;
    out.ampdist = Dy;
    out.f2_warped = f2_warped;
end
end

function [nll, g, aux] = geom_aux(geom, c)
[nll, g] = geom(c);
aux = [];
end

function SSE = warp_sse(c, Phi, rep, q1, q2, var1)
[~, ~, SSE] = bayes_warp_geom(c, Phi, rep, q1, q2, var1);
end

function out = statsFun(mat)
out = quantile(mat.',[0.025,.975]);
end

% GP prior N(0, s * (corr(L) + 0.1 I)) stored through its Cholesky factor
function gp = gp_prior(s, L, Dmat, Rcorr)
if nargin < 4
    Rcorr = chol(exp2corr2(L, Dmat) + 0.1*eye(size(Dmat,1)));
end
gp.Rcorr = Rcorr;
gp.s = s;
gp.logdet = 2*sum(log(diag(Rcorr))) + size(Dmat,1)*log(s);
end

function out = quad_form(R, f)
z = R.' \ f(:);
out = z.' * z;
end

function R = chol_jitter(K)
jit = 0;
n = size(K,1);
[R, p] = chol(K);
while p > 0
    jit = max(2*jit, 1e-10*mean(diag(K)));
    [R, p] = chol(K + jit*eye(n));
end
end

% Metropolis update of a latent function f (f1 or f2) given the data y, its
% GP prior and the SRVF alignment likelihood. f | y, prior is Gaussian with
% precision K^{-1} + I/sigmaf, so the target is that density times the
% alignment likelihood exp(-SSE/(2 sigma)).
function [f_curr, q_curr, SSE_curr, accept] = update_f(f_curr, q_curr, SSE_curr, y, sse_fun, gp, Rprop, t01, sigma, sigmaf)
f_prop = f_curr + randn(1, length(f_curr)) * Rprop;
q_prop = f_to_srvf(f_prop', t01, 0);
SSE_prop = sse_fun(q_prop);
lp_curr = -SSE_curr/(2*sigma) - (quad_form(gp.Rcorr, f_curr)/gp.s + sum((y - f_curr).^2)/sigmaf)/2;
lp_prop = -SSE_prop/(2*sigma) - (quad_form(gp.Rcorr, f_prop)/gp.s + sum((y - f_prop).^2)/sigmaf)/2;
accept = log(rand) < (lp_prop - lp_curr);
if accept
    f_curr = f_prop;
    q_curr = q_prop;
    SSE_curr = SSE_prop;
end
end

% Metropolis-Hastings update of a GP length scale with a truncated normal
% random-walk proposal (flat prior on L > 0)
function [L_curr, gp, accept] = update_length(f, gp, s, L_curr, propsd, Dmat)
L_prop = random(truncate(makedist('Normal', L_curr, propsd), 0, inf));
gp_prop = gp_prior(s, L_prop, Dmat);
lp_curr = -gp.logdet/2 - quad_form(gp.Rcorr, f)/(2*s);
lp_prop = -gp_prop.logdet/2 - quad_form(gp_prop.Rcorr, f)/(2*s);
% Hastings correction for the truncation: q(L'|L) = phi(.)/Phi(L/propsd)
log_hastings = log(normcdf(L_curr/propsd)) - log(normcdf(L_prop/propsd));
accept = log(rand) < (lp_prop - lp_curr + log_hastings);
if accept
    L_curr = L_prop;
    gp = gp_prop;
end
end

% Normalize gamma to [0,1]
function gam = norm_gam(gam)
gam = (gam-gam(1))./(gam(end)-gam(1));
end
