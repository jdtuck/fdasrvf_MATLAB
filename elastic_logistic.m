function out = elastic_logistic(f, y, t, option)
% ELASTIC_LOGISTIC Elastic Logistic Functional Regression
% -------------------------------------------------------------------------
% This function identifies a logistic regression model with 
% phase-variablity using elastic methods
%
% Usage:  out = elastic_logistic(f, y, t)
%         out = elastic_logistic(f, y, t, option)
% 
% Input:
% f (M,N): matrix defining N functions of M samples
% y : response vector of length N
% t : time vector of length M
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.B = []; % defines basis if empty uses bspline
% option.df = 20; % degress of freedom
% option.sparam = 25; % number of times to run filter
% option.max_itr = 20; % maximum number of iterations
%
% Output:
% structure with fields:
% alpha: intercept
% beta: regression function
% fn: aligned functions
% qn: aligned srvfs
% gamma: warping functions
% q: original srvfs
% B: basis Matrix used
% b: coefficient vector
% Loss: logistic loss
% type: model type

addpath(genpath('DP'))

if nargin < 4
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.B = [];
    option.df = 20;
    option.max_itr = 20;
end

if option.parallel == 1
    if isempty(gcp('nocreate'))
        % prompt user for number threads to use
        nThreads = input('Enter number of threads to use: ');
        if nThreads > 1
            parpool(nThreads);
        elseif nThreads > 12 % check if the maximum allowable number of threads is exceeded
            while (nThreads > 12) % wait until user figures it out
                fprintf('Maximum number of threads allowed is 12\n Enter a number between 1 and 12\n');
                nThreads = input('Enter number of threads to use: ');
            end
            if nThreads > 1
                parpool(nThreads);
            end
        end
    end
end

%% Initialize

a = size(t,1);
if (a ~=1)
    t = t';
end
binsize = mean(diff(t));
[M, N] = size(f);

if option.smooth == 1
    f = smooth_data(f, option.sparam);
end

% create B-spline basis
if isempty(option.B)
    B = create_basismatrix(t, option.df, 4);
else
    B = option.B;
end
Nb = size(B,2);

q = f_to_srvf(f,t);

gamma = repmat(linspace(0,1,M)',1,N);

%% Main Loop
itr = 1;
LL = zeros(1,option.max_itr);
while itr <= option.max_itr
    fprintf('Iteration: %d\n', itr);
    % align data
    fn = zeros(M,N);
    qn = zeros(M,N);
    for k = 1:N
        fn(:,k) = interp1(t, f(:,k), (t(end)-t(1)).*gamma(:,k) + t(1))';
        qn(:,k) = gradient(fn(:,k), binsize)./sqrt(abs(gradient(fn(:,k), binsize))+eps);
    end
    
    Phi = ones(N, Nb+1);
    for ii = 1:N
        for jj = 2:Nb+1
            Phi(ii,jj) = trapz(t, qn(:,ii) .* B(:,jj-1));
        end
    end
    
    % find alpha and beta using bfgs
    options.Method = 'lbfgs';
    options.Display = 'off';
    b0 = zeros(Nb+1, 1);
    b = minFunc(@logit_optim,b0,options,Phi,y);
    
    alpha = b(1);
    beta = B * b(2:Nb+1);
    
    % compute the lostic loss
    LL(itr) = logit_loss(b, Phi, y);
        
    % find gamma
    gamma_new = zeros(M,N);
    if option.parallel == 1
        parfor ii=1:N
            gamma_new(:,ii) = logistic_warp(beta, t, q(:,ii), y(ii));
        end
    else
        for ii=1:N
            gamma_new(:,ii) = logistic_warp(beta, t, q(:,ii), y(ii));
        end
    end
    
    if norm(gamma-gamma_new) < 1e-5
        break
    else
        gamma = gamma_new;
    end
    
    itr = itr + 1;
end
gamma = gamma_new;

out.alpha = alpha;
out.beta = beta;
out.fn = fn;
out.qn = qn;
out.gamma = gamma;
out.q = q;
out.B = B;
out.b = b(2:end);
out.SSE = LL(1:itr-1);
out.type = 'logistic';

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end

end

%% Helper Functions

function out = phi(t)
% calculates logisitc function, returns 1/(1+exp(-t))
idx = t > 0;
out = zeros(size(t));
out(idx) = 1./(1+exp(-t(idx)));
exp_t = exp(t(~idx));
out(~idx) = exp_t ./ (1+exp_t);
end

function out = logit_loss(b, X, y)
% logistic loss function, returns Sum{-log(phi(t))}
z = X * b;
yz = y.*z;
idx = yz > 0;
out = zeros(size(yz));
out(idx) = log(1+exp(-1.*yz(idx)));
out(~idx) = (-1.*yz(~idx) + log(1+exp(yz(~idx))));
out = sum(out);
end

function grad = logit_gradient(b, X, y)
% calculates gradient of the logistic loss
z = X * b;
z = phi(y.*z);
z0 = (z-1).*y;
grad = X.' * z0;
end

function Hs = logit_hessian(s, b, X, y)
% calculates hessian of the logistic loss
z = X * b;
z = phi(y.*z);
d = z.*(1-z);
wa = d.*(X*s);
Hs = X.' * wa;
end

function [nll, g] = logit_optim(b, X, y)
% function for call to optimizer
nll = logit_loss(b, X, y);
if nargout > 1
    g = logit_gradient(b, X, y);
end
end

function gamma = logistic_warp(beta, t, q, y)
% calculates optimal warping for functional logistic regression
q = q.';
if y == 1
    gamma = optimum_reparam(beta',q,t,0);
elseif y == -1
    gamma = optimum_reparam(-1.*beta',q,t,0);
end
end

