function out = pairwise_align_bayes(f1i, f2i, time, mcmcopts)
% PAIRWISE_ALIGN Align two functions using Bayesian method
% -------------------------------------------------------------------------
% This function aligns two functions using Bayesian framework. It will align
% f2 to f1. It is based on mapping warping functions to a hypersphere, and a
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
% Usage:  out = pairwise_align_bayes(f1i, f2i, time)
%         out = pairwise_align_bayes(f1i, f2i, time, mcmcopts)
%
% Input:
% f1i: vector defining M samples of function 1
% f2i: vector defining M samples of function 2
% time: time vector of length M
%
% default mcmc options
% mcmcopts.iter = 2e4; % number of iterations
% mcmcopts.burnin = min(5e3,mcmcopts.iter/2); % number of burnin
% mcmcopts.alpha0 = 0.1; # inverse gamma prior parameters
% mcmcopts.beta0 = 0.1; # inverse gamma prior parameters
% tmp.betas = [0.5,0.5,0.005,0.0001]; %pczn paramaters
% tmp.probs = [0.1,0.1,0.7,0.1];
% mcmcopts.zpcn = tmp;
% mcmcopts.propvar = 1; % proposal vairance
% mcmcopts.initcoef = repelem(0, 20).'; % init coef
% mcmcopts.npoints = 200; % number of sample interpolation points
% mcmcopts.extrainfo = true; % return extra info about mcmc
%
% Output:
% structure containing
% out.f2_warped: aligned f2
% out.gamma: warping function
% out.g_coef: final g_coef
% out.psi: final psi
% out.sigma1: final sigma
%
% if extrainfo
% out.accept: accept of psi samples
% out.betas_ind
% out.logl: log likelihood
% out.gamma_mat: posterior gammas
% out.gamma_stats: posterior gamma stats
% out.xdist: phase distance posterior
% out.ydist: amplitude distance posterior

if nargin < 4
    mcmcopts.iter = 2e4;
    mcmcopts.burnin = min(5e3,mcmcopts.iter/2);
    mcmcopts.alpha0 = 0.1;
    mcmcopts.beta0 = 0.1;
    tmp.betas = [0.5,0.5,0.005,0.0001];
    tmp.probs = [0.1,0.1,0.7,0.1];
    mcmcopts.zpcn = tmp;
    mcmcopts.propvar = 1;
    mcmcopts.initcoef = repelem(0, 20).';
    mcmcopts.npoints = 200;
    mcmcopts.extrainfo = true;
end

if (~iscolumn(f1i))
    f1i = f1i.';
end
if (~iscolumn(f2i))
    f2i = f2i.';
end
if (~iscolumn(time))
    time = time.';
end
if (length(f1i) ~= length(f2i))
    error('Length of f1 and f2 must be equal')
end
if (length(f1i) ~= length(time))
    error('Length of f1 and time must be equal')
end
if (length(mcmcopts.zpcn.betas) ~= length(mcmcopts.zpcn.probs))
    error('In zpcn, betas must equal length of probs')
end
if (mod(length(mcmcopts.initcoef), 2) ~= 0)
    error('Length of mcmcopts.initcoef must be even')
end

% Number of sig figs to report in gamma_mat
SIG_GAM = 13;
iter = mcmcopts.iter;

% for now, back to struct format of Yi's software
f1.x = time;
f1.y = f1i;
f2.x = time;
f2.y = f2i;

% normalize timet to [0,1]
% ([a,b] - a) / (b-a) = [0,1]
f1.x = (f1.x-min(f1.x))./(max(f1.x)-min(f1.x));
f2.x = (f2.x-min(f2.x))./(max(f2.x)-min(f2.x));

% parameter settings
pw_sim_global_burnin = mcmcopts.burnin;
valid_index = pw_sim_global_burnin:iter;
pw_sim_global_Mg = length(mcmcopts.initcoef)/2;
g_coef_ini = mcmcopts.initcoef;
numSimPoints = mcmcopts.npoints;
pw_sim_global_domain_par = linspace(0,1,numSimPoints).';
g_basis = basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mg, 1);
sigma1_ini = 1;
zpcn = mcmcopts.zpcn;
pw_sim_global_sigma_g = mcmcopts.propvar;

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

% srsf transformation
q1 = f_Q(f1);
q2 = f_Q(f2);

tmp = f_exp1(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis, false));
if (min(tmp.y)<0)
    error("Invalid initial value of g")
end

% result objects
result.g_coef = zeros(iter,length(g_coef_ini));
result.sigma1 = zeros(1,iter);
result.logl = zeros(1,iter);
result.SSE = zeros(1,iter);
result.accept = zeros(1,iter);
result.accept_betas = zeros(1,iter);

% init
g_coef_curr = g_coef_ini;
sigma1_curr = sigma1_ini;
SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false),q1,q2);
logl_curr = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false),q1,q2,sigma1_ini^2,SSE_curr);

result.g_coef(1,:) = g_coef_ini;
result.sigma1(1) = sigma1_ini;
result.SSE(1) = SSE_curr;
result.logl(1) = logl_curr;

% update the chain for iter-1 times
obj = ProgressBar(iter-1, ...
    'Title', 'Running MCMC' ...
    );
for m = 2:iter
    % update g
    [g_coef_curr, ~, SSE_curr, accept, zpcnInd] = f_updateg_pw(g_coef_curr, g_basis, sigma1_curr^2, q1, q2, SSE_curr, @propose_g_coef);

    % update sigma1
    newshape = length(q1.x)/2 + mcmcopts.alpha0;
    newscale = 1/2 * SSE_curr + mcmcopts.beta0;
    sigma1_curr = sqrt(1/gamrnd(newshape,1/newscale));
    logl_curr = f_logl_pw(f_basistofunction(g_basis.x, 0, g_coef_curr, g_basis, false), q1, q2, sigma1_curr^2, SSE_curr);

    % save update to results
    result.g_coef(m,:) = g_coef_curr;
    result.sigma1(m) = sigma1_curr;
    result.SSE(m) = SSE_curr;
    if (mcmcopts.extrainfo)
        result.logl(m) = logl_curr;
        result.accept(m) = accept;
        result.accept_betas(m) = zpcnInd;
    end
    obj.step([], [], []);
end
obj.release();

% calculate posterior mean of psi
pw_sim_est_psi_matrix = zeros(length(pw_sim_global_domain_par), length(valid_index));
for k = 1:length(valid_index)
    g_temp = f_basistofunction(pw_sim_global_domain_par, 0, result.g_coef(valid_index(k),:).', g_basis, false);
    psi_temp = f_exp1(g_temp);
    pw_sim_est_psi_matrix(:,k) = psi_temp.y;
end

result_posterior_psi_simDomain = f_psimean(pw_sim_global_domain_par, pw_sim_est_psi_matrix);

% resample to same number of points as the input f1 and f2
result_i = interp1(result_posterior_psi_simDomain.x, result_posterior_psi_simDomain.y, f1.x, 'linear', 'extrap');
result_posterior_psi.x=f1.x;
result_posterior_psi.y=result_i;

% transform posterior mean of psi to gamma
result_posterior_gamma = f_phiinv(result_posterior_psi);
gam0 = result_posterior_gamma.y;
result_posterior_gamma.y = norm_gam(gam0);

% warped f2
f2_warped = warp_f_gamma(f2.y, result_posterior_gamma.y, result_posterior_gamma.x);

if (mcmcopts.extrainfo)
    % matrix of posterior draws from gamma
    gamma_mat = zeros(length(q1.x),size(pw_sim_est_psi_matrix,2));
    one_v = ones(1,size(pw_sim_est_psi_matrix,1));
    Dx = zeros(1,size(pw_sim_est_psi_matrix,2));
    Dy = Dx;
    for ii = 1:size(pw_sim_est_psi_matrix,2)
        result_i = interp1(result_posterior_psi_simDomain.x, pw_sim_est_psi_matrix(:,ii), f1.x, 'linear', 'extrap');
        result_i2.y=result_i;
        result_i2.x=f1.x;
        tmp = f_phiinv(result_i2);
        gamma_mat(:,ii) = round(norm_gam(tmp.y),SIG_GAM);
        v = inv_exp_map(one_v, pw_sim_est_psi_matrix(:,ii));
        Dx(ii) = sqrt(trapz(pw_sim_global_domain_par, v.^2));
        q2warp = warp_q_gamma(q2.y, gamma_mat(:,ii), q2.x).';
        Dy(ii) = sqrt(trapz(q2.x,(q1.y-q2warp).^2)); 
    end
    gamma_stats = statsFun(gamma_mat);
end

% return object
out.f2_warped = f2_warped;
out.gamma = result_posterior_gamma.y;
out.g_coef = result.g_coef;
out.psi = result_posterior_psi.y;
out.sigma1 = result.sigma1;

if (mcmcopts.extrainfo)
    out.accept = result.accept(2:end);
    out.betas_ind = result.accept_betas(2:end);
    out.logl = result.logl;
    out.gamma_mat = gamma_mat;
    out.gamma_stats = gamma_stats.';
    out.xdist = Dx;
    out.ydist = Dy;
end
end

function out = statsFun(vec)
a = quantile(vec.',0.025);
b = quantile(vec.',0.975);
out = [a;b];
end

function out = f_exp1(g)
out.x = g.x;
out.y = bcalcY(f_L2norm(g), g.y);
end

% function for calculating the next MCMC sample given current state
function [g_coef, logl, SSE, accept, zpcnInd] = f_updateg_pw(g_coef_curr,g_basis,var1_curr,q1,q2,SSE_curr,propose_g_coef)
g_coef_prop = propose_g_coef(g_coef_curr);

tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
while (min(tst.y)<0)
    g_coef_prop = propose_g_coef(g_coef_curr);
    tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
end

if (SSE_curr == 0)
    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis, false), q1, q2);
end

SSE_prop = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis,false), q1, q2);

logl_curr = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q1, q2, var1_curr, SSE_curr);

logl_prop = f_logl_pw(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis,false), q1, q2, var1_curr, SSE_prop);

ratio = min(1, exp(logl_prop-logl_curr));

u = rand;
if (u <= ratio)
    g_coef = g_coef_prop.prop;
    logl = logl_prop;
    SSE = SSE_prop;
    accept = true;
    zpcnInd = g_coef_prop.ind;
end

if (u > ratio)
    g_coef = g_coef_curr;
    logl = logl_curr;
    SSE = SSE_curr;
    accept = false;
    zpcnInd = g_coef_prop.ind;
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
function out = f_SSEg_pw(g, q1, q2)
obs_domain = q1.x;
exp1g_temp = f_predictfunction(f_exp1(g), obs_domain, 0);
pt = [0; bcuL2norm2(obs_domain, exp1g_temp.y)];
tmp = f_predictfunction(q2, pt, 0);
vec = (q1.y - tmp.y .* exp1g_temp.y).^2;
out = sum(vec);
end

function out = f_logl_pw(g, q1, q2, var1, SSEg)
if (SSEg == 0)
    SSEg = f_SSEg_pw(g, q1, q2);
end
n = length(q1.y);
out = n * log(1/sqrt(2*pi)) - n * log(sqrt(var1)) - SSEg / (2 * var1);
end

%##########################################################################
% calculate Q(f), Qinv(q)
% f,q: function in the form of list$x, list$y
% fini: f(0)
% returns Q(f),Qinv(q), function in the form of struct.x, struct.y
%##########################################################################
function out = f_Q(f)
d = f_predictfunction(f, f.x, 1);
out.x = f.x;
out.y = sign(d.y) .* sqrt(abs(d.y));
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
