function out = elastic_regression(f, y, t, lambda, option)
% ELASTIC_REGRESSION Elastic Functional Regression
% -------------------------------------------------------------------------
% This function identifies a regression model with phase-variablity using
% elastic methods
%
% Usage:  out = elastic_regression(f, y, t)
%         out = elastic_regression(f, y, t, lambda)
%         out = elastic_regression(f, y, t, lambda, option)
%
% input:
% f (M,N): matrix defining N functions of M samples
% y : response vector of length N
% t : time vector of length M
% lambda: regularization parameter
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 1; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.B = []; % defines basis if empty uses bspline
% option.df = 20; % degress of freedom
% option.sparam = 25; % number of times to run filter
% option.max_itr = 20; % maximum number of iterations
%
% output:
% structure with fields:
% alpha: intercept
% beta: regression function
% fn: aligned functions
% qn: aligned srvfs
% gamma: warping functions
% q: original srvfs
% B: basis Matrix used
% b: coefficient vector
% SSE: sum of squared errors
% type: model type

addpath(genpath('DP'))

if nargin < 4
    lambda = 0;
    option.parallel = 0;
    option.closepool = 1;
    option.smooth = 0;
    option.sparam = 25;
    option.B = [];
    option.df = 20;
    option.max_itr = 20;
elseif nargin < 5
    option.parallel = 0;
    option.closepool = 1;
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

%% Parameters

fprintf('\n lambda = %5.1f \n', lambda);

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

% second derivative for regularization
Bdiff = zeros(M, Nb);
for ii = 1:Nb
    Bdiff(:,ii) = gradient(gradient(B(:,ii), binsize),binsize);
end

q = f_to_srvf(f,t);

gamma = repmat(linspace(0,1,M)',1,N);

itr = 1;
SSE = zeros(1,option.max_itr);

while itr <= option.max_itr
    fprintf('Iteration: %d\n', itr);
    % align data
    fn = zeros(M,N);
    qn = zeros(M,N);
    for k = 1:N
        fn(:,k) = interp1(t, f(:,k), (t(end)-t(1)).*gamma(:,k) + t(1))';
        qn(:,k) = gradient(fn(:,k), binsize)./sqrt(abs(gradient(fn(:,k), binsize))+eps);
    end
    
    % OLS using basis
    Phi = ones(N, Nb+1);
    for ii = 1:N
        for jj = 2:Nb+1
            Phi(ii,jj) = trapz(t, qn(:,ii) .* B(:,jj-1));
        end
    end
    
    R = zeros(Nb+1, Nb+1);
    for ii = 2:Nb+1
        for jj = 2:Nb+1
            R(ii,jj) = trapz(t, Bdiff(:,ii-1).*Bdiff(:,ii-1));
        end
    end
    
    xx = Phi.' * Phi;
    inv_xx = xx + lambda * R;
    xy = Phi.' * y;
    b = inv_xx\xy;
    
    alpha = b(1);
    beta = B * b(2:Nb+1);
    
    % compute the SSE
    int_X = zeros(N,1);
    for ii = 1:N
        int_X(ii) = trapz(t, qn(:,ii).*beta);
    end
    
    SSE(itr) = sum((y-alpha-int_X).^2);
    
    % find gamma
    gamma_new = zeros(M,N);
    if option.parallel == 1
        parfor ii=1:N
            gamma_new(:,ii) = regression_warp(beta, t, q(:,ii), y(ii), alpha);
        end
    else
        for ii=1:N
            gamma_new(:,ii) = regression_warp(beta, t, q(:,ii), y(ii), alpha);
        end
    end

    if norm(gamma-gamma_new) < 1e-5
        break
    else
        gamma = gamma_new;
    end
    
    itr = itr + 1;
end

% last step with centering of gamma
gamI = SqrtMeanInverse(gamma_new.');
gamI_dev = gradient(gamI, 1/(M-1));
beta = interp1(t, beta, (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
for k = 1:N
    qn(:,k) = interp1(t, qn(:,k), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
    fn(:,k) = interp1(t, fn(:,k), (t(end)-t(1)).*gamI + t(1))';
    gamma(:, k) = interp1(t, gamma_new(:, k), (t(end)-t(1)).*gamI + t(1));
end

out.alpha = alpha;
out.beta = beta;
out.fn = fn;
out.qn = qn;
out.gamma = gamma;
out.q = q;
out.B = B;
out.b = b(2:end);
out.SSE = SSE(1:itr-1);
out.type = 'linear';

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end
end

%% Helper Functions
function gamma = zero_crossing(Y, q, bt, t, y_max,y_min, gmax, gmin)
% finds zero-crossing of optimal gamma, gam = s*gmax + (1-s)*gmin
% from elastic regression model
max_itr = 100;
M = length(t);
a = zeros(1, max_itr);
a(1) = 1;
f = zeros(1, max_itr);
f(1) = y_max - Y;
f(2) = y_min - Y;
mrp = f(1);
mrn = f(2);
mrp_ind = 1;  % most recent positive index
mrn_ind = 2;  % most recent negative index

for ii = 3:max_itr
    x1 = a(mrp_ind);
    x2 = a(mrn_ind);
    y1 = mrp;
    y2 = mrn;
    a(ii) = (x1*y2 - x2*y1) / (y2-y1);
    
    gam_m = a(ii) * gmax + (1 - a(ii)) * gmin;
    gamdev = gradient(gam_m, 1/(M-1));
    qtmp = interp1(t, q, (t(end)-t(1)).*gam_m + t(1)).*sqrt(gamdev); 
    f(ii) = trapz(t, qtmp.*bt) - Y;
    
    if abs(f(ii)) < 1e-5
        break
    elseif f(ii) > 0
        mrp = f(ii);
        mrp_ind = ii;
    else
        mrn = f(ii);
        mrn_ind = ii;
    end
end

gamma = a(ii) * gmax + (1 - a(ii)) * gmin;

end

function gamma_new = regression_warp(beta, t, q, y, alpha)
% calculates optimal warping for function linear regression
M = length(t);
beta = beta';
q = q';
gam_M = optimum_reparam(beta,q,t,0);
gamdev = gradient(gam_M, 1/(M-1));
qM = interp1(t, q, (t(end)-t(1)).*gam_M + t(1)).*sqrt(gamdev);
y_M = trapz(t, qM.*beta);

gam_m = optimum_reparam(-1.*beta,q,t,0);
gamdev = gradient(gam_m, 1/(M-1));
qm = interp1(t, q, (t(end)-t(1)).*gam_m + t(1)).*sqrt(gamdev);
y_m = trapz(t, qm.*beta);

if y > alpha + y_M
    gamma_new = gam_M;
elseif y < alpha + y_m
    gamma_new = gam_m;
else
    gamma_new = zero_crossing(y-alpha, q, beta, t, y_M, y_m, gam_M, gam_m);
end

end
