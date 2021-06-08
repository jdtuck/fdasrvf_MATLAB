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
% mcmcopts.h = .01; % stepsize
% mcmcopts.L = 4;  % number of leapfrog steps
% mcmcopts.f1propvar = 0.0001;
% mcmcopts.f2propvar = 0.0001; % proposal variance for f2
% mcmcopts.L1propvar = 0.3; %controls L1 (in f1 cov) acceptance rate
% mcmcopts.L2propvar = 0.3;  %controls L2 (in f2 cov) acceptance rate
% mcmcopts.npoints = 200;
% mcmcopts.vpriorvar = 1; % proposal variance for theta
% mcmcopts.initcoef = repelem(0, 20).'; % init coef
% mcmcopts.nbasis = 10; %number of basis of gradient g
% mcmcopts.basis = 'fourier';  %'fourier' 'legendre'
% mcmcopts.extrainfo = true;
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
%
% J. D. Tucker, L. Shand, and K. Chowdhary. “Multimodal Bayesian Registration of Noisy Functions using 
%       Hamiltonian Monte Carlo”, Computational Statistics and Data Analysis, accepted, 2021.
%

if nargin < 4
    mcmcopts.iter = 1e4;
    mcmcopts.burnin = min(5e3,mcmcopts.iter/2); %must be >1
    mcmcopts.thin = 50;
    mcmcopts.nchains = 1;
    mcmcopts.alpha0 = 0.1;
    mcmcopts.beta0 = 0.1;
    mcmcopts.alpha = 1;
    mcmcopts.beta = 1;
    mcmcopts.h = .01; % stepsize
    mcmcopts.L = 4;  % number of leapfrog steps
    mcmcopts.f1propvar = 0.0001;
    mcmcopts.f2propvar = 0.0001; % proposal variance for f2
    mcmcopts.L1propvar = 0.3; %controls L1 (in f1 cov) acceptance rate
    mcmcopts.L2propvar = 0.3;  %controls L2 (in f2 cov) acceptance rate
    %mcmcopts.npoints = 200;
    mcmcopts.vpriorvar = 1; % proposal variance for theta
    mcmcopts.initcoef = repelem(0, 20).'; % init coef
    mcmcopts.nbasis = 10; %number of basis of gradient g
    mcmcopts.basis = 'fourier';  %'fourier' 'legendre'
    mcmcopts.sampfreq = 5; %controls smoothness of fitrgp
    mcmcopts.seed = randi(floor(intmax/10));
    mcmcopts.extrainfo = true;
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
if (mod(mcmcopts.nbasis, 2) ~= 0)
    error('Length of mcmcopts.nbasis must be even')
end

% setup random start points for more than 1 chain
random_starts = zeros(length(mcmcopts.initcoef),mcmcopts.nchains);
if mcmcopts.nchains > 1
    rng(mcmcopts.seed)
    for i = 1:mcmcopts.nchains
        randcoef = -1 + (2)*rand(length(mcmcopts.initcoef),1); %.1*(2*rand(length(mcmcopts.initcoef),1) - 1.0);
        random_starts(:,i) = randcoef;
    end
end

isparallel = true;
try
    tmp = ver('parallel');
catch
    isparallel = false;
end

if mcmcopts.nchains == 1
    isparallel = false;
end

if isparallel
    nThreads = min(feature('numcores'),mcmcopts.nchains);
    if isempty(gcp('nocreate'))
        parpool(nThreads);
    end
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
            psi1 = sqrt(gradient(out.gamma(:,i),binsize));
            psi2 = sqrt(gradient(out.gamma(:,j),binsize));
            Dx(i,j) = acos(trapz(time1,psi1.*psi2));
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
        [~,gam_mu,~,~] = SqrtMean(tmp.');
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

% for now, back to struct format of Yi's software
y1.x = time;
y1.y = y1i;
y2.x = time;
y2.y = y2i;
T = length(time);

% normalize timet to [0,1]
% ([a,b] - a) / (b-a) = [0,1]
y1.x = (y1.x-min(y1.x))./(max(y1.x)-min(y1.x));
y2.x = (y2.x-min(y2.x))./(max(y2.x)-min(y2.x));

% warping function parameter settings
pw_sim_global_burnin = mcmcopts.burnin;
valid_index = pw_sim_global_burnin:iter;
ncoef = length(mcmcopts.initcoef); %number of basis functions for warping function v
nbasis = mcmcopts.nbasis; %number of basis functions for d_basis in HMC
pw_sim_global_Mv = ncoef/2;
numSimPoints = T;
pw_sim_global_domain_par = linspace(0,1,numSimPoints).';
d_basis = basis_fourierd(pw_sim_global_domain_par, nbasis);
if strcmpi(mcmcopts.basis,'fourier')
    v_basis = basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mv, 1);
elseif strcmpi(mcmcopts.basis,'legendre')
    v_basis = basis_legendre(pw_sim_global_domain_par, pw_sim_global_Mv, 1);
else
    error('Incorrect Basis Specified')
end
sigma_ini = 1;
v_priorvar = mcmcopts.vpriorvar;
v_coef_ini = mcmcopts.initcoef;
D = pdist(time');
Dmat = squareform(D);
C = v_priorvar ./ repelem(1:pw_sim_global_Mv,2);
cholC = chol(diag(C));
h = mcmcopts.h;
L = mcmcopts.L;

    function result = propose_v_coef(v_coef_curr)
        v_coef_new = normrnd(v_coef_curr, C', ncoef, 1);
        result = v_coef_new;
    end

% f1,f2 prior, propoposal params
sigma1_ini = 0.01;
sigma2_ini = 0.01;
f1_propvar = mcmcopts.f1propvar;
f2_propvar = mcmcopts.f2propvar;
%model1 = fitrgp(time', y1.y');
%model2 = fitrgp(time', y2.y');
model1 = fitrgp(time(:,1:mcmcopts.sampfreq:end)', y1.y(1:mcmcopts.sampfreq:end)');
model2 = fitrgp(time(:,1:mcmcopts.sampfreq:end)', y2.y(1:mcmcopts.sampfreq:end)');

s1_ini = model1.KernelInformation.KernelParameters(2)^2; %variance of f1, s_1^2, shouldn we initialize randomly?
s2_ini = model2.KernelInformation.KernelParameters(2)^2; %variance of f1, s_2^2
L1_propvar = mcmcopts.L1propvar;
L2_propvar = mcmcopts.L2propvar;
L1_ini = model1.KernelInformation.KernelParameters(1);
L2_ini = model2.KernelInformation.KernelParameters(1);

K_f1_corr = exp2corr2(L1_ini,Dmat)+0.1 * eye(length(y1.y));
K_f1 = s1_ini * K_f1_corr;
K_f1 = inv(K_f1);
K_f2_corr = exp2corr2(L2_ini,Dmat)+ 0.1 * eye(length(y2.y));
K_f2 = s2_ini * K_f2_corr;
K_f2 = inv(K_f2);
K_f1prop= exp2corr(f1_propvar,L1_ini,Dmat);
K_f2prop= exp2corr(f2_propvar,L2_ini,Dmat);

% result objects
result.v_coef = zeros(iter,length(v_coef_ini));
result.sigma = zeros(1,iter);
result.sigma1 = zeros(1,iter);
result.sigma2 = zeros(1,iter);
result.f1 = zeros(iter,length(time));
result.f2_warped = result.f1;
result.f2 = zeros(iter,length(time));
result.s1 = zeros(1,iter);
result.s2 = zeros(1,iter);
result.L1 = zeros(1,iter);
result.L2 = zeros(1,iter);

result.logl = zeros(1,iter);
result.SSE = zeros(1,iter);
result.theta_accept = zeros(1,iter);
result.f2_accept = zeros(1,iter);
result.f1_accept = zeros(1,iter);
result.L2_accept = zeros(1,iter);
result.L1_accept = zeros(1,iter);

% init
v_coef_curr = v_coef_ini;
v_curr =  f_basistofunction(v_basis.x,0 , v_coef_curr, v_basis, false);
sigma_curr = sigma_ini;
sigma1_curr = sigma1_ini;
sigma2_curr = sigma2_ini;
L1_curr = L1_ini;
L2_curr = L2_ini;

%f1_curr = resubPredict(model1)'; %these serve as proposal means so must be smooth
%f2_curr = resubPredict(model2)';
f1_curr = predict(model1,time').';
f2_curr = predict(model2,time').';


% SRVF transformation setup
q1_curr.x = y1.x;
q2_curr.x = y2.x;
q1_curr.y = f_to_srvf(f1_curr.',time,0);
q2_curr.y = f_to_srvf(f2_curr.',time,0);

SSE_curr = f_SSEv_pw(v_curr,q1_curr,q2_curr);
logl_curr = f_vpostlogl_pw(v_curr, q1_curr, q2_curr, sigma_curr, SSE_curr);

result.v_coef(1,:) = v_coef_ini;
result.f1(1,:) = f1_curr;
result.f2(1,:) = f2_curr;
result.f2_warped(1,:) = f2_curr;
result.sigma(1) = sigma_ini;
result.sigma1(1) = sigma1_ini;
result.sigma2(1) = sigma2_ini;
result.s1(1) = s1_ini;
result.s2(1) = s2_ini;
result.L1(1) = L1_ini;
result.L2(1) = L2_ini;
result.SSE(1) = SSE_curr;
result.SSEprop(1) = SSE_curr;
result.logl(1) = logl_curr;

if mcmcopts.nchains == 1
    % update the chain for iter-1 times
    obj = ProgressBar(iter-1, ...
        'Title', 'Running MCMC' ...
        );
end
[M,N] = size(f1_curr);
n = M*N;

[nll, g, SSE_curr] = f_dlogl_pw(v_coef_curr, v_basis, d_basis, sigma_curr, q1_curr, q2_curr);
for m = 2:iter
    
    % update f1
    [f1_curr,q1_curr, f1_accept] = f_updatef1_pw(f1_curr,q1_curr, y1,q2_curr,v_coef_curr, v_basis,SSE_curr,K_f1,K_f1prop,sigma_curr,sqrt(sigma1_curr));
    
    % update f2
    [f2_curr,q2_curr, f2_accept] = f_updatef2_pw(f2_curr,q2_curr, y2,q1_curr,v_coef_curr, v_basis,SSE_curr,K_f2,K_f2prop,sigma_curr,sqrt(sigma2_curr));
    
    % update v
    [v_coef_curr, nll, g, SSE_curr, theta_accept] = f_updatev_pw(v_coef_curr, v_basis, sqrt(sigma_curr), q1_curr, q2_curr,nll, g,SSE_curr,@propose_v_coef,d_basis,cholC,h,L);
    
    % update sigma^2
    newshape = length(q1_curr.x)/2 + mcmcopts.alpha;
    newscale = 1/2 * SSE_curr + mcmcopts.beta;
    sigma_curr = 1/gamrnd(newshape,1/newscale);  % @todo should we paramatize as precision?
    
    % update sigma1^2
    newshape = n/2 + mcmcopts.alpha0;
    newscale = sum((y1.y - f1_curr).^2)/2 + mcmcopts.beta0;
    sigma1_curr = 1/gamrnd(newshape,1/newscale);
    
    % update sigma2^2
    newshape = n/2 + mcmcopts.alpha0;
    newscale = sum((y2.y - f2_curr).^2)/2 + mcmcopts.beta0;
    sigma2_curr = 1/gamrnd(newshape,1/newscale);
    
    %update hyerparameters
    
    % update s1^2
    newshape = n/2 + mcmcopts.alpha0;
    newscale = (f1_curr/K_f1_corr * f1_curr.')/2 + mcmcopts.beta0;
    s1_curr = 1/gamrnd(newshape,1/newscale);
    
    % update s2^2
    newshape = n/2 + mcmcopts.alpha0;
    newscale = (f2_curr/K_f2_corr * f2_curr.')/2 + mcmcopts.beta0;
    s2_curr = 1/gamrnd(newshape,1/newscale);
    
    % update L1
    [L1_curr, L1_accept] = f_updatephi_pw(f1_curr,K_f1,s1_curr, L1_curr, L1_propvar, Dmat);
    
    % update L2
    [L2_curr, L2_accept] = f_updatephi_pw(f2_curr,K_f2,s2_curr, L2_curr, L2_propvar, Dmat);
    
    
    K_f1_corr = exp2corr2(L1_curr,Dmat) + 0.1 * eye(length(y1.y));
    K_f1 = s1_curr * K_f1_corr;
    K_f1 = inv(K_f1);
    K_f2_corr = exp2corr2(L2_curr,Dmat)+ 0.1 * eye(length(y2.y));
    K_f2 = s2_curr * K_f2_corr;
    K_f2 = inv(K_f2);
    
    v_curr =  f_basistofunction(v_basis.x,0 , v_coef_curr, v_basis, false);
    logl_curr = f_vpostlogl_pw(v_curr, q1_curr, q2_curr, sigma_curr, SSE_curr);
    
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
    result.logl(m) = logl_curr;
    if (mcmcopts.extrainfo)
        result.theta_accept(m) = theta_accept;
        result.f1_accept(m) = f1_accept;
        result.f2_accept(m) = f2_accept;
        result.L1_accept(m) = L1_accept;
        result.L2_accept(m) = L2_accept;
    end
    if mcmcopts.nchains == 1
        obj.step([], [], []);
    end
end
if mcmcopts.nchains == 1
    obj.release();
end

% calculate posterior mean of psi
pw_sim_est_psi_matrix = zeros(length(pw_sim_global_domain_par), length(valid_index));
for k = 1:length(valid_index)
    v_temp = f_basistofunction(pw_sim_global_domain_par, 0, result.v_coef(valid_index(k),:).', v_basis, false);
    psi_temp = f_exp1(v_temp);
    pw_sim_est_psi_matrix(:,k) = psi_temp.y;
end

result_posterior_psi_simDomain = f_psimean(pw_sim_global_domain_par, pw_sim_est_psi_matrix);

% resample to same number of points as the input f1 and f2
result_i = interp1(result_posterior_psi_simDomain.x, result_posterior_psi_simDomain.y, y1.x, 'linear', 'extrap');
result_posterior_psi.x=y1.x.';
result_posterior_psi.y=result_i.';

% transform posterior mean of psi to gamma
result_posterior_gamma = f_phiinv(result_posterior_psi); %fix this
gam0 = result_posterior_gamma.y;
result_posterior_gamma.y = norm_gam(gam0);

if (mcmcopts.extrainfo)
    % matrix of posterior draws from gamma
    gamma_mat = zeros(length(y1.x),size(pw_sim_est_psi_matrix,2));
    one_v = ones(1,size(pw_sim_est_psi_matrix,1));
    Dx = zeros(1,size(pw_sim_est_psi_matrix,2));
    Dy = Dx;
    for ii = 1:size(pw_sim_est_psi_matrix,2)
        result_i = interp1(psi_temp.x, pw_sim_est_psi_matrix(:,ii), y1.x, 'linear', 'extrap');
        result_i2.y=result_i.';
        result_i2.x=y1.x.';
        tmp = f_phiinv(result_i2);
        gamma_mat(:,ii) = round(norm_gam(tmp.y),SIG_GAM);
        v = inv_exp_map(one_v, pw_sim_est_psi_matrix(:,ii));
        Dx(ii) = sqrt(trapz(pw_sim_global_domain_par, v.^2));
        q2warp = warp_q_gamma(q2_curr.y, gamma_mat(:,ii), y2.x).';
        Dy(ii) = sqrt(trapz(y2.x,(q1_curr.y-q2warp).^2));
    end
    gamma_stats = statsFun(gamma_mat);
end


% return object
out.f1 = result.f1(valid_index,:);
out.f2 = result.f2(valid_index,:);
out.gamma = result_posterior_gamma.y;
out.v_coef = result.v_coef(valid_index,:);
out.psi = result_posterior_psi.y;
out.sigma = result.sigma(valid_index);
out.sigma1 = result.sigma1(valid_index);
out.sigma2 = result.sigma2(valid_index);
out.s1 = result.s1(valid_index);
out.s2 = result.s2(valid_index);
out.L1 = result.L1(valid_index);
out.L2 = result.L2(valid_index);
out.f2_warped_mu = warp_f_gamma(mean(out.f2,1), result_posterior_gamma.y, result_posterior_gamma.x);


if (mcmcopts.extrainfo)
    out.theta_accept = result.theta_accept(2:end);
    out.f1_accept = result.f1_accept(2:end);
    out.f2_accept = result.f2_accept(2:end);
    out.L1_accept = result.L1_accept(2:end);
    out.L2_accept = result.L2_accept(2:end);
    out.gamma_mat = gamma_mat;
    out.gamma_stats = gamma_stats.';
    out.SSE = result.SSE;
    out.logl = result.logl;
    out.phasedist = Dx;
    out.ampdist = Dy;
    
    f2_warped = zeros(length(valid_index),length(result_posterior_gamma.x));
    for ii = 1:length(valid_index)
        f2_warped(ii,:) = warp_f_gamma(out.f2(ii,:)', gamma_mat(:,ii), result_posterior_gamma.x);
    end
    out.f2_warped = f2_warped;
end
end

function out = statsFun(mat)
out = quantile(mat.',[0.025,.975]);
end

function out = f_exp1(v)
out.x = v.x;
out.y = bcalcY(f_L2norm(v), v.y);
end

% updating theta/v
function [v_coef_curr,nll, g, SSE_curr, theta_accept] = f_updatev_pw(v_coef_curr,v_basis,sigma_curr,q1,q2,nll_cur,g_cur,SSE_curr,propose_v_coef,d_basis,cholC,h,L) %sigma_curr=variance
v_coef_prop = propose_v_coef(v_coef_curr);

%qfunc = f_basistofunction(v_basis.x,0 , v_coef_prop.prop', v_basis, false);
q = v_coef_prop; %proposed coefficients for warping function
D = length(q);
rth = sqrt(h); % make the scale comparable to MALA
g_cur = -1*g_cur; %g_cur should be current set of basis coefficients

% sample velocity
v = randn(D,1);
v = cholC'*v;

% natural gradient
halfng = cholC*g_cur;
ng = cholC*halfng;

% accumulate the power of force
pow = rth/2*(g_cur'*v);

% calculate current energy
E_cur = nll_cur - h/8*(halfng'*halfng);

randL = ceil(rand*L);

% Alternate full sth for position and velocity
for l = 1:randL
    % Make a half step for velocity
    v = v + rth/2*ng;
    
    % Make a full step for position
    rot = (q+1i*v).*exp(-1i*rth);
    q = real(rot); v = imag(rot);
    
    % update geometry
    qvec.y=q;
    qvec.x = q1.x';
    [nll, g, SSEv] = f_dlogl_pw(qvec.y,v_basis, d_basis, sigma_curr, q1, q2);
    halfng = cholC*g;
    ng = cholC'*halfng;
    
    % Make a half step for velocity
    v = v + rth/2*ng;
    
    % accumulate the power of force
    if l~=randL
        pow = pow + rth*(g'*v);
    end
end

% accumulate the power of force
pow = pow + rth/2*(g'*v);

% Evaluate energy at start and end of trajectory
E_prp = nll - h/8*(halfng'*halfng);

% Accept or reject the state at end of trajectory, returning either
% the position at the end of the trajectory or the initial position
logRatio = - E_prp + E_cur - pow;

if (isfinite(logRatio) && (log(rand) < min([0,logRatio])))
    v_coef_curr = q;
    theta_accept = true;
    SSE_curr = SSEv;
else
    nll = nll_cur; g = g_cur;
    theta_accept = false;
end
end

%update f2
function [f2_curr,q2_curr, f2_accept] = f_updatef2_pw(f2_curr,q2_curr,y2,q1,v_coef_curr, v_basis, SSE_curr,K_f2,K_f2prop,sigma_curr,sigma2_curr)
time = y2.x;
v = f_basistofunction(v_basis.x,0,v_coef_curr,v_basis,false);

f2_prop = mvnrnd(f2_curr, K_f2prop);
q2_prop.x = y2.x;
q2_prop.y = f_to_srvf(f2_prop',time,0);

SSE_prop = f_SSEv_pw(v, q1, q2_prop);

postlog_curr = f_f2postlogl_pw(f2_curr,y2,SSE_curr,K_f2,sigma_curr,sigma2_curr);
postlog_prop = f_f2postlogl_pw(f2_prop,y2,SSE_prop,K_f2,sigma_curr,sigma2_curr);

ratio = min(1, exp(postlog_prop-postlog_curr));

u = rand;
if (u <= ratio)
    f2_curr = f2_prop;
    q2_curr = q2_prop;
    f2_accept = true;
else
    f2_accept = false;
end
end

%update f1
function [f1_curr,q1_curr, f1_accept] = f_updatef1_pw(f1_curr,q1_curr,y1,q2,v_coef_curr, v_basis, SSE_curr,K_f1,K_f1prop,sigma_curr,sigma1_curr)
time = y1.x;
v = f_basistofunction(v_basis.x,0,v_coef_curr,v_basis,false);

f1_prop = mvnrnd(f1_curr, K_f1prop);
q1_prop.x = y1.x;
q1_prop.y = f_to_srvf(f1_prop',time,0);

SSE_prop = f_SSEv_pw(v, q1_prop, q2);

postlog_curr = f_f1postlogl_pw(f1_curr,y1,SSE_curr,K_f1,sigma_curr,sigma1_curr);
postlog_prop = f_f1postlogl_pw(f1_prop,y1,SSE_prop,K_f1,sigma_curr,sigma1_curr);

ratio = min(1, exp(postlog_prop-postlog_curr));

u = rand;
if (u <= ratio)
    f1_curr = f1_prop;
    q1_curr = q1_prop;
    f1_accept = true;
else
    f1_accept = false;
end
end

%update L1, L2
function [L1_curr, L1_accept] = f_updatephi_pw(f1_curr,K_f1,s1_curr,L1_curr, L1_propvar, Dmat)
pd = makedist('Normal',L1_curr, L1_propvar);
pd_trunc = truncate(pd, 0,inf);
L1_prop = random(pd_trunc);

K_f1_tmp = s1_curr * (exp2corr2(L1_prop,Dmat) +0.1 * eye(length(f1_curr)));

SSEf_curr = (f1_curr * K_f1 * f1_curr.')/2;
SSEf_prop = (f1_curr/K_f1_tmp * f1_curr.')/2;

%K_f1 is inverse but K_f1_tmp is not
postlog_prop = -log(det(K_f1_tmp))/2 - SSEf_prop - log(1-normcdf(-L1_curr/L1_propvar));
postlog_curr = log(det(K_f1))/2 - SSEf_curr - log(1-normcdf(-L1_prop/L1_propvar));
ratio = min(1, exp(postlog_prop-postlog_curr));

u = rand;
if (u <= ratio)
    L1_curr = L1_prop;
    L1_accept = true;
else
    L1_accept = false;
end
end

%##########################################################################
% For pairwise registration, evaluate the neg log likelihood of v, given q1 and
% q2
% v is given in the form of struct.x and struct.y
% q1 and q2 are vectors
% var1: model variance
% SSEv: if not provided, than re-calculate
% returns a numeric value which is nlogl(v|q1,q2), see Eq 10 of JCGS
%##########################################################################
% SSEv: sum of sq error= sum over ti of { q1(ti)-{q2, v}(ti) }^2

function out = f_warp_pw(v, q1, q2)
obs_domain = q1.x;
exp1v_temp = f_predictfunction(f_exp1(v), obs_domain, 0);
pt = [0; bcuL2norm2(obs_domain.', exp1v_temp.y.')];
tmp = f_predictfunction(q2, pt, 0);
out = tmp.y.' .* exp1v_temp.y;
end

function out = f_SSEv_pw(v, q1, q2)
q2_gamma = f_warp_pw(v,q1,q2);
vec = (q1.y - q2_gamma.').^2;
out = sum(vec);
end

%Posterior Log Likelihood of v
function [out, SSEv] = f_vpostlogl_pw(v, q1, q2, var, SSEv)
if (SSEv == 0)
    SSEv = f_SSEv_pw(v, q1, q2);
end
n = length(q1.y);
out = - n * log(sqrt(var)) - SSEv / (2 * var);
end

%##########################################################################
% Related Functions to posterior of f1 and f2 and respective covaraince params
%##########################################################################
% SSEv: sum of sq error= sum over ti of { q1(ti)-{q2, v}(ti) }^2

%Posterior Log Likelihood of f1
function out = f_f1postlogl_pw(f1,y1,SSEv,K_f1,sigma_curr,sigma1_curr)
n = length(y1.x);
iSig_f1 = K_f1+ eye(n) .* n/sigma1_curr;
f1_mean = iSig_f1 \ y1.y.' .* n/sigma1_curr;
SSEf1 = (f1-f1_mean.') * iSig_f1 * (f1.'-f1_mean);
out = -SSEv/(2 * sigma_curr) - SSEf1/2;
end

%Posterior Log Likelihood of f2
function out = f_f2postlogl_pw(f2,y2,SSEv,K_f2,sigma_curr,sigma2_curr)
n = length(y2.x);
iSig_f2 = K_f2 + eye(n) .* n/sigma2_curr;
f2_mean = iSig_f2 \ y2.y.' .* n/sigma2_curr;
SSEf2 = (f2-f2_mean.') * iSig_f2 * (f2.'-f2_mean);
out = -SSEv/(2 * sigma_curr) - SSEf2/2;
end

%##########################################################################
% Extrapolate a function given by a discrete vector
% f: function in the form of list$x, list$y
% at: t values for which f(t) is returned
% deriv: can calculate derivative
% method: smoothing method: 'linear' (default), 'cubic'
% returns: $y==f(at), $x==at
%##########################################################################
function out = f_predictfunction(f, at, deriv)
if (deriv == 0)
    result = interp1(f.x,f.y,at,'linear','extrap');
    out.x = at;
    out.y = result;
end

if (deriv == 1)
    fmod = interp1(f.x,f.y,at,'linear','extrap');
    diffy1 = [0; diff(fmod)];
    diffy2 = [diff(fmod); 0];
    diffx1 = [0; diff(at)];
    diffx2 = [diff(at); 0];
    
    
    out.x = at;
    out.y = (diffy2 + diffy1) ./ (diffx2 + diffx1);
end
end

%##########################################################################
% calculate L2 norm of a function, using trapezoid rule for integration
% f:function in the form of list$x, list$y
% returns ||f||, a numeric value
%##########################################################################
function out = f_L2norm(f)
out = border_l2norm(f.x,f.y);
end

%##########################################################################
% Different basis functions b_i()
% f.domain: grid on which b_i() is to be returned
% numBasis: numeric value, number of basis functions used
% (note: #basis = #coef/2 for Fourier basis)
% fourier.p: period of the Fourier basis used
% returns a struct:
%     matrix: with nrow=length(t) and ncol=numBasis (or numBasis*2 for
%     Fourier)
%     x: f.domain
%##########################################################################
function out = basis_fourier(f_domain, numBasis, fourier_p)
result = zeros(length(f_domain), 2*numBasis);
for i = 1:(2*numBasis)
    j = ceil(i/2);
    if (mod(i,2) == 1)
        result(:,i) = sqrt(2) * sin(2*j*pi*f_domain./fourier_p);
    end
    if (mod(i,2) == 0)
        result(:,i) = sqrt(2) * cos(2*j*pi*f_domain./fourier_p);
    end
    out.x = f_domain;
    out.matrix = result;
end
end

%##########################################################################
% Legendre (non normalized) basis functions b_i()
% f.domain: grid on which b_i() is to be returned
% numBasis: numeric value, number of basis functions used
% (note: #basis = #coef/2 for legendre basis)
% fourier.p: period of the Fourier basis used
% returns a struct:
%     matrix: with nrow=length(t) and ncol=numBasis (or numBasis*2 for
%     Fourier)
%     x: f.domain
%##########################################################################
function out = basis_legendre(f_domain, numBasis, fourier_p)
result = zeros(length(f_domain), 2*numBasis);
for i = 1:2*numBasis
    f_domain_scaled = 2*(f_domain./fourier_p) - 1;
    tmp = legendre(i,f_domain_scaled);
    result(:,i) = tmp(1,:);
end
out.x = f_domain;
out.matrix = result;
end

%##########################################################################
% Given the coefficients of basis functions, returns the actual function on
% a grid
% f.domain: numeric vector, grid of the actual function to return
% coefconst: leading constant term
% coef: numeric vector, coefficients of the basis functions
%       Note: if #coef < #basis functions, only the first %coef basis
%             functions will be used
% basis: in the form of list$x, list$matrix
% plotf: if true, show a plot of the function generated
% returns the generated function in the form of struct.x=f.domain, struct.y
%##########################################################################
function result = f_basistofunction(f_domain, coefconst, coef, basis, plotf)
if (size(basis.matrix,2) < length(coef))
    error('In f_basistofunction, #coeffients exceeds #basis functions.')
end
result.x = basis.x;
result.y = basis.matrix(:,(1:length(coef))) * coef + coefconst;

result = f_predictfunction(result, f_domain, 0);
if (plotf)
    plot(result.x,result.y)
end

end

%##########################################################################
% Calculate Karcher mean/median with Alg1/Alg2 in (Kurtek,2014)
% x: vector of length = length(domain of the psi's)
% y: M columns, each of length = length(x)
% e1, e2: small positive constants
% method: 'ext' = extrinsic, 'int' = intrinsic
% returns posterier mean/median of M psi's (of form .x, .y)
%##########################################################################
function out = f_psimean(x, y)
rmy = mean(y,2);
tmp.x = x;
tmp.y = rmy;
result = rmy / f_L2norm(tmp);
out.x = x;
out.y = result;
end

function out = f_phiinv(psi)
f_domain = psi.x;
result = [0; bcuL2norm2(f_domain, psi.y)];
out.x = f_domain;
out.y = result;
end

% Normalize gamma to [0,1]
function gam = norm_gam(gam)
gam = (gam-gam(1))./(gam(end)-gam(1));
end

%##########################################################################
% For pairwise registration, evaluate the derivative loglikelihood of
% gradient g for coefficients of basis rep of warping function, given q1 and q2
% q1, q2 are all given in the form of struct.x and struct.y
% g is a vector
% var1: model variance
% returns a numeric value which is dglogl(g|q1,q2)
%##########################################################################
function [nll, g_coef, SSEv] = f_dlogl_pw(v_coef, v_basis, d_basis, sigma_curr, q1, q2)
vec = f_basistofunction(v_basis.x,0 , v_coef, v_basis, false);
psi = f_exp1(vec);
N = length(q1.x);
obs_domain = linspace(0,1,N)';
binsize = mean(diff(obs_domain));
gamma = f_phiinv(psi);
q2_warp = warp_q_gamma(q2.y, gamma.y,obs_domain,1)';
q2_warp_grad = gradient(q2_warp,binsize);

basismat = d_basis.matrix;

g = zeros(N,1);
for i = 1:size(basismat,2)
    ubar = cumtrapz(obs_domain,basismat(:,i));
    integrand = (q1.y-q2_warp).*(-2.*q2_warp_grad.*ubar-q2_warp.*basismat(:,i));
    tmp = 1/sigma_curr * trapz(obs_domain,integrand);
    g = g + tmp.*basismat(:,i);
end

[out, SSEv] = f_vpostlogl_pw(vec, q1, q2, sigma_curr, 0);

nll = -1*out;
g_coef = v_basis.matrix' * g;
end
