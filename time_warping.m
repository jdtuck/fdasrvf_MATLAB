function out = time_warping(f,t,lambda,option)
% Group-wise function alignment
%
% This function aligns a collection of functions using the elastic square-root
% slope (srsf) framework.
% 
% input:
% f (M,N): matrix defining N functions of M samples
% t : time vector of length M
% lambda: regularization parameter
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting
% option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
% option.w = 0.0; % BFGS weight
% option.MaxItr = 20;  % maximum iterations
%
% output structure containing
% fn: aligned functions
% qn: aligned srvfs
% q0: original srvfs
% fmean: function mean
% mqn: mean srvf
% gam: warping functions
% psi: srvf of gam
% stats: structure of statistics of alignment
addpath(genpath('DP'))

if nargin < 3
    lambda = 0;
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
elseif nargin < 4
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
end

% time warping on a set of functions
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

% check dimension of time
a = size(t,1);
if a == 1
    t = t';
end

binsize = mean(diff(t));
[M, N] = size(f);
f0 = f;

if option.smooth == 1
    f = smooth_data(f, option.sparam);
end

if option.showplot == 1
    figure(1); clf;
    plot(t, f, 'linewidth', 1);
    title('Original data', 'fontsize', 16);
    pause(0.1);
end

%% Compute the q-function of the plot
q = f_to_srvf(f,t);

%% Set initial using the original f space
fprintf('\nInitializing...\n');
mnq = mean(q,2);
dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
[~, min_ind] = min(dqq);
mq = q(:,min_ind);
mf = f(:,min_ind);

gam = zeros(N,size(q,1));
if option.parallel == 1
    parfor k = 1:N
        q_c = q(:,k,1); mq_c = mq;
        gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                   mf(1), f(1,k,1));
    end
else
    for k = 1:N
        q_c = q(:,k,1); mq_c = mq;
        gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                   mf(1), f(1,k,1));
    end
end
gamI = SqrtMeanInverse(gam);
mf = warp_f_gamma(mf,gamI,t);
mq = f_to_srvf(mf,t);

%% Compute Mean
fprintf('Computing Karcher mean of %d functions in SRVF space...\n',N);
ds = inf;
MaxItr = option.MaxItr;
qun = zeros(1,MaxItr);
f_temp = zeros(length(t),N);
q_temp = zeros(length(t),N);
for r = 1:MaxItr
    fprintf('updating step: r=%d\n', r);
    if r == MaxItr
        fprintf('maximal number of iterations is reached. \n');
    end

    % Matching Step
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    gam = zeros(N,size(q,1));
    gam_dev = zeros(N,size(q,1));
    if option.parallel == 1
        parfor k = 1:N
            q_c = q(:,k,1); mq_c = mq(:,r);
            gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                       mf(1,r), f(1,k,1));
            gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
            f_temp(:,k) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
            q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
        end
    else
        for k = 1:N
            q_c = q(:,k,1); mq_c = mq(:,r);
            gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                       mf(1,r), f(1,k,1));
            gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
            f_temp(:,k) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
            q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
        end
    end
    q(:,:,r+1) = q_temp;
    f(:,:,r+1) = f_temp;

    ds(r+1) = sum(simps(t, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
        lambda*sum(simps(t, (1-sqrt(gam_dev')).^2));

    % Minimization Step
    % compute the mean of the matched function
    mq(:,r+1) = mean(q(:,:,r+1),2);
    mf(:,r+1) = mean(f(:,:,r+1),2);

    qun(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
    if qun(r) < 1e-2 || r >= MaxItr
        break;
    end
end

% last step with centering of gam
r = r+1;
if option.parallel == 1
    parfor k = 1:N
        q_c = q(:,k,1); mq_c = mq(:,r);
        gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                   mf(1,r), f(1,k,1));
    end
else
    for k = 1:N
        q_c = q(:,k,1); mq_c = mq(:,r);
        gam(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                                   mf(1,r), f(1,k,1));
    end
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
mq(:,r+1) = interp1(t, mq(:,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
for k = 1:N
    q(:,k,r+1) = interp1(t, q(:,k,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
    f(:,k,r+1) = interp1(t, f(:,k,r), (t(end)-t(1)).*gamI + t(1))';
    gam(k,:) = interp1(t, gam(k,:), (t(end)-t(1)).*gamI + t(1));
end

%% Aligned data & stats
fn = f(:,:,r+1);
qn = q(:,:,r+1);
q0 = q(:,:,1);
mean_f0 = mean(f0, 2);
std_f0 = std(f0, 0, 2);
mean_fn = mean(fn, 2);
std_fn = std(fn, 0, 2);
mqn = mq(:,r+1);
fmean = mean(f0(1,:))+cumtrapz(t,mqn.*abs(mqn));

fgam = zeros(M,N);
for ii = 1:N
    fgam(:,ii) = interp1(t, fmean, (t(end)-t(1)).*gam(ii,:) + t(1));
end
var_fgam = var(fgam,[],2);

stats.orig_var = trapz(t,std_f0.^2);
stats.amp_var = trapz(t,std_fn.^2);
stats.phase_var = trapz(t,var_fgam);

gam = gam.';
[~,fy] = gradient(gam,binsize,binsize);
psi = sqrt(fy+eps);

if option.showplot == 1
    figure(2); clf;
    plot((0:M-1)/(M-1), gam, 'linewidth', 1);
    axis square;
    title('Warping functions', 'fontsize', 16);

    figure(3); clf;
    plot(t, fn, 'LineWidth',1);
    title(['Warped data, \lambda = ' num2str(lambda)], 'fontsize', 16);

    figure(4); clf;
    plot(t, mean_f0, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_f0+std_f0, 'r-', 'linewidth', 1);
    plot(t, mean_f0-std_f0, 'g-', 'linewidth', 1);
    title('Original data: Mean \pm STD', 'fontsize', 16);

    figure(5); clf;
    plot(t, mean_fn, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_fn+std_fn, 'r-', 'linewidth', 1);
    plot(t, mean_fn-std_fn, 'g-', 'linewidth', 1);
    title(['Warped data, \lambda = ' num2str(lambda) ': Mean \pm STD'], 'fontsize', 16);

    figure(6); clf;
    plot(t, fmean, 'g','LineWidth',1);
    title(['f_{mean}, \lambda = ' num2str(lambda)], 'fontsize', 16);
end

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end

out.f0 = f;
out.time = t;
out.fn = fn;
out.qn = qn;
out.q0 = q0;
out.fmean = fmean;
out.mqn = mqn;
out.gam = gam;
out.psi = psi;
out.stats = stats;
out.qun = qun(1:r);
out.lambda = lambda;
out.method = option.method;
out.gamI = gamI;
out.rsamps = false;
