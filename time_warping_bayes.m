function obj = time_warping_bayes(f, time, mcmcopts)
% TIME_WARPING_BAYES Align multiple functions using Bayesian method
% -------------------------------------------------------------------------
% This function aligns a collection of functions using the elastic square-root
% slope (srsf) framework. This method uses a Bayesian approach.
% It is based on mapping warping functions to a hypersphere, and a
% subsequent exponential mapping to a tangent space. In the tangent space,
% the Z-mixture pCN algorithm is used to explore both local and global
% structure in the posterior distribution.
%
% The Z-mixture pCN algorithm uses a mixture distribution for the proposal
% distribution, controlled by input parameter zpcn. The zpcn$betas must be
% between 0 and 1, and are the coefficients of the mixture components, with
% larger coefficients corresponding to larger shifts in parameter space. The
% zpcn$probs give the probability of each shift size.
%
% Usage:  out = time_warping_bayes(f, time)
%         out = time_warping_bayes(f, time, mcmcopts)
%
% Input:
% f: (M,N): matrix defining N functions of M samples
% time: time vector of length M
%
% default mcmc options
% mcmcopts.iter = 2000; % number of iterations
% mcmcopts.burnin = min(5e3,mcmcopts.iter/2); % number of burnin
% mcmcopts.alpha0 = 0.1; # inverse gamma prior parameters
% mcmcopts.beta0 = 0.1; # inverse gamma prior parameters
% tmp.betas = [0.5,0.5,0.005,0.0001]; %pczn paramaters
% tmp.probs = [0.1,0.1,0.7,0.1];
% mcmcopts.zpcn = tmp;
% mcmcopts.propvar = 1; % proposal variance
% mcmcopts.initcoef = zeros(20, size(f,2)) % init warping function coef
% mcmcopts.nbasis = 40; % number of basis elements for q
% mcmcopts.npoints = 200; % number of sample interpolation points
% mcmcopts.extrainfo = true; % return extra info about mcmc
%
% Output:
% fdawarp object
%
% if extrainfo object contains
% mcmc.accept: accept of q samples
% mcmc.betas_ind
% mcmc.gamma_mat: posterior gammas
% mcmc.gamma_stats: posterior gamma stats
[M, N] = size(f);

if nargin < 4
    mcmcopts.iter = 2000;
    mcmcopts.burnin = min(5e3,mcmcopts.iter/2);
    mcmcopts.alpha0 = 0.1;
    mcmcopts.beta0 = 0.1;
    tmp.betas = [0.5,0.5,0.005,0.0001];
    tmp.probs = [0.1,0.1,0.7,0.1];
    mcmcopts.zpcn = tmp;
    mcmcopts.propvar = 1;
    mcmcopts.initcoef = zeros(20, N);
    mcmcopts.nbasis = 40;
    mcmcopts.npoints = 200;
    mcmcopts.extrainfo = true;
end

if (~iscolumn(time))
    time = time.';
end

if (size(f,1) ~= length(time))
    error('Length of f1 and time must be equal')
end

if (length(mcmcopts.zpcn.betas) ~= length(mcmcopts.zpcn.probs))
    error('In zpcn, betas must equal length of probs')
end

if (mod(size(mcmcopts.initcoef), 1) ~= 0)
    error('Length of mcmcopts.initcoef must be even')
end

if (mod(mcmcopts.nbasis, 1) ~= 0)
    error('mcmcopts.nbasis must be even')
end

% Number of sig figs to report in gamma_mat
SIG_GAM = 13;
iter = mcmcopts.iter;

% normalize timet to [0,1]
% ([a,b] - a) / (b-a) = [0,1]
time = (time-min(time))./(max(time)-min(time));

% parameter settings
pw_sim_global_burnin = mcmcopts.burnin;
valid_index = pw_sim_global_burnin:iter;
pw_sim_global_Mg = size(mcmcopts.initcoef,1)/2;
g_coef_ini = mcmcopts.initcoef;
numSimPoints = mcmcopts.npoints;
pw_sim_global_domain_par = linspace(0,1,numSimPoints).';
g_basis = basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mg, 1);
pw_sim_global_Mq = mcmcopts.nbasis/2;
g_basis_q = basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mq, 1);
sigma1_ini = 1;
zpcn = mcmcopts.zpcn;
pw_sim_global_sigma_g = mcmcopts.propvar;
pw_sim_global_sigma_q = mcmcopts.propvar;

    function result = propose_g_coef(g_coef_curr)
        pCN_beta = zpcn.betas;
        pCN_prob = zpcn.probs;
        probm = [0, cumsum(pCN_prob)];
        z = rand;
        for i = 1:length(pCN_beta)
            if (z <= probm(i+1) && z > probm(i))
                g_coef_new = normrnd(0, pw_sim_global_sigma_g ./ repelem(1:pw_sim_global_Mg,2), 1, pw_sim_global_Mg * 2);
                result.prop = sqrt(1-pCN_beta(i)^2) * g_coef_curr + pCN_beta(i) * g_coef_new.';
                result.ind = i;
            end
        end
    end

    function result = propose_q_star(q_star_curr)
        pCN_beta = zpcn.betas;
        pCN_prob = zpcn.probs;
        probm = [0, cumsum(pCN_prob)];
        z = rand;
        for i = 1:length(pCN_beta)
            if (z <= probm(i+1) && z > probm(i))
                q_star_new = normrnd(0, pw_sim_global_sigma_q ./ repelem(1:pw_sim_global_Mq,2), 1, pw_sim_global_Mq * 2);
                result.prop = sqrt(1-pCN_beta(i)^2) * q_star_curr + pCN_beta(i) * q_star_new.';
                result.ind = i;
            end
        end
    end

% srsf transformation
f0 = zeros(length(g_basis.x),N);
for ii = 1:N
    f0(:,ii) = interp1(time,f(:,ii),g_basis.x,'linear');
end
qo = f_to_srvf(f0, g_basis.x);


for ii = 1:N
    tmp = f_exp1(f_basistofunction(g_basis.x,0,g_coef_ini(:,ii),g_basis, false));
    if (min(tmp.y)<0)
        error("Invalid initial value of g for sample %d", ii)
    end
end

% result objects
result.g_coef = zeros(iter,size(g_coef_ini,1),size(g_coef_ini,2));
result.q_star = zeros(iter,length(g_basis.x));
result.sigma1 = zeros(1,iter);
result.logl = zeros(1,iter);
result.SSE = zeros(iter,N);
result.accept = zeros(1,iter);
result.accept_betas = zeros(1,iter);

% init
g_coef_curr = g_coef_ini;
sigma1_curr = sigma1_ini;
mnq = mean(qo,2);
dqq = sqrt(sum((qo - mnq*ones(1,N)).^2,1));
[~, min_ind] = min(dqq);
q_star_curr = qo(:,min_ind);
q.y = qo;
q.x = g_basis.x;
q_star_ini = q_star_curr;
SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false), ...
    q,q_star_ini);
logl_curr = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false), ...
    q,q_star_ini,sigma1_ini^2,SSE_curr);

result.g_coef(1,:,:) = g_coef_ini;
result.sigma1(1) = sigma1_ini;
result.q_star(1,:) = q_star_curr;
result.SSE(1,:) = SSE_curr;
result.logl(1) = logl_curr;

q_star_coef_curr = g_basis_q.matrix\q_star_curr;
clear q_star_curr
q_star_curr.x = g_basis.x;
q_star_curr.y = q_star_ini;

% update the chain for iter-1 times
barobj = ProgressBar(iter-1, ...
    'Title', 'Running MCMC' ...
    );
for m = 2:iter
    % update g
    for ii = 1:N
        q_tmp = q;
        q_tmp.y = q_tmp.y(:,ii);
        [g_coef_curr(:,ii), ~, ~] = f_updateg_pw(g_coef_curr(:,ii), g_basis, sigma1_curr^2, q_tmp, q_star_curr.y, @propose_g_coef);
    end
    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
        q,q_star_curr.y);
    % update q_star
    [q_star_coef_curr, accept, zpcnInd] = f_updateq_pw(g_coef_curr, g_basis_q, sigma1_curr^2, q, q_star_coef_curr, SSE_curr, @propose_q_star);
    q_star_curr = f_basistofunction(g_basis_q.x,0,q_star_coef_curr,g_basis_q,false);
    
    % center
    result_posterior_gamma = zeros(length(g_basis.x),N);
    for ii = 1:N
        g_temp = f_basistofunction(pw_sim_global_domain_par, 0, g_coef_curr(:,ii), g_basis, false);
        psi_temp = f_exp1(g_temp);
        % transform posterior mean of psi to gamma
        tmp_gamma = f_phiinv(psi_temp);
        gam0 = tmp_gamma.y;
        tmp_gamma.y = norm_gam(gam0);
        result_posterior_gamma(:,ii) = norm_gam(gam0);
    end
    gamI_o = SqrtMeanInverse(result_posterior_gamma.');
    q_star_curr.y = warp_q_gamma(q_star_curr.y,gamI_o,q_star_curr.x);
    for k = 1:N
        result_posterior_gamma(:,k) = interp1(q_star_curr.x, result_posterior_gamma(:,k), (time(end)-time(1)).*gamI_o + time(1));
        gamma.x = g_basis.x;
        gamma.y = result_posterior_gamma(:,k);
        psi = f_phi(gamma);
        [~,yy] = f_exp1inv(psi);
        g_coef_curr(:,k) = g_basis.matrix\yy;
    end
    
    % update sigma1
    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
        q,q_star_curr.y);
    newshape = N*length(time)/2 + mcmcopts.alpha0;
    newscale = 1/2 * sum(SSE_curr) + mcmcopts.beta0;
    sigma1_curr = sqrt(1/gamrnd(newshape,1/newscale));
    logl_curr = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
        q,q_star_curr.y,sigma1_curr^2,SSE_curr);
    
    % save update to results
    result.g_coef(m,:,:) = g_coef_curr;
    result.q_star(m,:) = q_star_curr.y;
    result.sigma1(m) = sigma1_curr;
    result.SSE(m,:) = SSE_curr;
    if (mcmcopts.extrainfo)
        result.logl(m) = logl_curr;
        result.accept(m) = accept;
        result.accept_betas(m) = zpcnInd;
    end
    barobj.step([], [], []);
end
barobj.release();

% calculate posterior mean of psi
pw_sim_est_psi_matrix = zeros(length(pw_sim_global_domain_par), length(valid_index),N);
for k = 1:length(valid_index)
    for ii = 1:N
        g_temp = f_basistofunction(pw_sim_global_domain_par, 0, result.g_coef(valid_index(k),:,ii).', g_basis, false);
        psi_temp = f_exp1(g_temp);
        pw_sim_est_psi_matrix(:,k,ii) = psi_temp.y;
    end
end

result_posterior_psi_simDomain = cell(N,1);
result_posterior_psi = zeros(M,N);
result_posterior_gamma = zeros(M,N);
f_warped = zeros(M,N);
for ii = 1:N
    result_posterior_psi_simDomain{ii} = f_psimean(pw_sim_global_domain_par, pw_sim_est_psi_matrix(:,:,ii));
    % resample to same number of points as the input f1 and f2
    result_i = interp1(result_posterior_psi_simDomain{ii}.x, result_posterior_psi_simDomain{ii}.y, time, 'linear', 'extrap');
    tmp_psi.x=time;
    tmp_psi.y=result_i;
    result_posterior_psi(:,ii) = tmp_psi.y;
    % transform posterior mean of psi to gamma
    tmp_gamma = f_phiinv(tmp_psi);
    gam0 = tmp_gamma.y;
    tmp_gamma.y = norm_gam(gam0);
    result_posterior_gamma(:,ii) = norm_gam(gam0);
    f_warped(:,ii) = warp_f_gamma(f(:,ii), result_posterior_gamma(:,ii), tmp_gamma.x);
end

if (mcmcopts.extrainfo)
    % matrix of posterior draws from gamma
    gamma_mat = zeros(length(time),size(pw_sim_est_psi_matrix,2),N);
    gamma_stats = zeros(2,size(pw_sim_est_psi_matrix,2),N);
    for jj = 1:N
        for ii = 1:size(pw_sim_est_psi_matrix,2)
            result_i = interp1(result_posterior_psi_simDomain{jj}.x, pw_sim_est_psi_matrix(:,ii,jj), time, 'linear', 'extrap');
            result_i2.y=result_i;
            result_i2.x=time;
            tmp = f_phiinv(result_i2);
            gamma_mat(:,ii,jj) = round(norm_gam(tmp.y),SIG_GAM);
            gamma_stats(:,ii,jj) = statsFun(gamma_mat(:,ii,jj));
        end
    end
end

% return object
obj.fn = f_warped;
obj.time = time;
obj.gam = result_posterior_gamma;
obj.psi = result_posterior_psi;
obj.qn = f_to_srvf(f_warped,time);
obj.q0 = f_to_srvf(f,time);
obj.mqn = mean(obj.qn,2);
obj.fmean = mean(f(1,:))+cumtrapz(obj.time,obj.mqn.*abs(obj.mqn));
std_f0 = std(f, 0, 2);
std_fn = std(obj.fn, 0, 2);
fgam = zeros(M,N);
for ii = 1:N
    fgam(:,ii) = warp_f_gamma(obj.fmean,obj.gam(:,ii),obj.time);
end
var_fgam = var(fgam,[],2);

obj.stats.orig_var = trapz(obj.time,std_f0.^2);
obj.stats.amp_var = trapz(obj.time,std_fn.^2);
obj.stats.phase_var = trapz(obj.time,var_fgam);

obj.qun = result.logl;
obj.method = "bayesian";
obj.gamI = interp1(g_basis.x, gamI_o, time, 'linear', 'extrap');
obj.rsamps = false;
obj.type = 'mean';

obj.mcmc.g_coef = result.g_coef;
obj.mcmc.sigma1 = result.sigma1;

if (mcmcopts.extrainfo)
    obj.mcmc.accept = result.accept(2:end);
    obj.mcmc.betas_ind = result.accept_betas(2:end);
    obj.mcmc.gamma_mat = gamma_mat;
    obj.mcmc.gamma_stats = gamma_stats;
end

end

function out = statsFun(vec)
a = quantile(vec,0.025);
b = quantile(vec,0.975);
out = [a,b];
end

function out = f_exp1(g)
out.x = g.x;
out.y = bcalcY(f_L2norm(g), g.y);
end

function [x,yy] = f_exp1inv(psi)
x = psi.x;
y = psi.y;

[x,ia,~] = unique(x);
[x, ia1] = sort(x);
y = y(ia);
y = y(ia1);
inner = round(trapzCpp(x,y), 10);

if ((inner < 1.001) && (inner >= 1))
    inner = 1;
end
if ((inner <= -1) && (inner > -1.001))
    inner = -1;
end
if ((inner < (-1)) || (inner > 1))
    fprintf("exp1inv: can't calculate the acos of: %f\n", inner);
end

theta = acos(inner);
yy = theta / sin(theta) .* (y - repelem(inner,length(y))');

if (theta==0)
    yy = zeros(1,length(x));
end

end

% function for calculating the next MCMC sample given current state
function [g_coef, accept, zpcnInd] = f_updateg_pw(g_coef_curr,g_basis,var1_curr,q1,q2,propose_g_coef)
g_coef_prop = propose_g_coef(g_coef_curr);

tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
while (min(tst.y)<0)
    g_coef_prop = propose_g_coef(g_coef_curr);
    tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
end

SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis, false), q1, q2);

SSE_prop = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis,false), q1, q2);

logl_ratio = -(SSE_prop - SSE_curr)/(2*var1_curr);

ratio = min(1, exp(logl_ratio));

u = rand;
if (u <= ratio)
    g_coef = g_coef_prop.prop;
    accept = true;
    zpcnInd = g_coef_prop.ind;
end

if (u > ratio)
    g_coef = g_coef_curr;
    accept = false;
    zpcnInd = g_coef_prop.ind;
end
end

% function for calculating the next MCMC sample given current state
function [q_star, accept, zpcnInd] = f_updateq_pw(g_coef_curr,g_basis,var1_curr,q,q_star_coef_curr,SSE_curr,propose_g_coef)
q_star_coef_prop = propose_g_coef(q_star_coef_curr);
ind = q_star_coef_prop.ind;
q_star_prop = f_basistofunction(g_basis.x,0,q_star_coef_prop.prop,g_basis, false);
q_star_curr = f_basistofunction(g_basis.x,0,q_star_coef_curr,g_basis, false);

if (SSE_curr == 0)
    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis, false), q, q_star_curr);
end

SSE_prop = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_prop.y);

logl_curr = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_curr.y, var1_curr, SSE_curr);

logl_prop = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_prop.y, var1_curr, SSE_prop);

ratio = min(1, exp(logl_prop-logl_curr));

u = rand;
if (u <= ratio)
    q_star = q_star_coef_prop.prop;
    accept = true;
    zpcnInd = ind;
end

if (u > ratio)
    q_star = q_star_coef_curr;
    accept = false;
    zpcnInd = ind;
end
end

%##########################################################################
% For pairwise registration, evaluate the loglikelihood of g, given q1 and
% q2
% g, q1, q2 are all given in the form of struct.x and struct.y
% var1: model variance
% SSEg: if not provided, than re-calculate
% returns a numeric value which is logl(g|q1,q2), see Eq 10 of JCGS
%##########################################################################
% SSEg: sum of sq error= sum over ti of { q1(ti)-{q2,g}(ti) }^2
% (Eq 11 of JCGS)
function out = f_SSEg_pw(g, q, q_star)
obs_domain = g.x;
out = zeros(1,size(q.y,2));
for ii = 1:size(q.y,2)
    g_tmp = g;
    g_tmp.y = g_tmp.y(:,ii);
    exp1g_temp = f_predictfunction(f_exp1(g_tmp), obs_domain, 0);
    pt = [0; bcuL2norm2(obs_domain, exp1g_temp.y)];
    q_tmp = q;
    q_tmp.y = q_tmp.y(:,ii);
    tmp = f_predictfunction(q_tmp, pt, 0);
    vec = (tmp.y .* exp1g_temp.y - q_star).^2;
    out(ii) = sum(vec);
end
end

function out = f_logl_pw(g, q, q_star, var1, SSEg)
if (SSEg == 0)
    SSEg = f_SSEg_pw(g, q, q_star);
end
n = length(q.x);
out = n * log(1/sqrt(2*pi)) + n * log(1/sqrt(var1)) - (SSEg ./ (2 * var1));
out = sum(out);
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
    out.y = zeros(length(at),size(f.y,2));
    out.x = at;
    for ii = 1:size(f.y,2)
        result = interp1(f.x,f.y(:,ii),at,'linear','extrap');
        out.y(:,ii) = result;
    end
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
if (size(basis.matrix,2) < size(coef,1))
    error('In f_basistofunction, #coeffients exceeds #basis functions.')
end
result.x = basis.x;
result.y = basis.matrix(:,(1:size(coef,1))) * coef + coefconst;

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

%##########################################################################
% calculate phi(gamma), phiinv(psi)
% gamma, psi: function in the form of struct.x, struct.y
% returns phi(gamma) or phiinv(psi), function in the form of struct.x,
% struct.y
%##########################################################################
function result = f_phi(gamma)
f_domain = gamma.x;
k = f_predictfunction(gamma, f_domain, 1);
k = k.y;
if (isempty(find(k < 0, 1)) ~= 0)
    idx = k < 0;
    k(idx) = 0;
end
result.x = f_domain;
result.y = sqrt(k);
if (f_L2norm(result) >= (1.01) || f_L2norm(result) <= (0.99))
    result.y = result.y / f_L2norm(result);
end
end
% the function returned by phi(gamma) = psi is always positive and has L2norm 1

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
