function out = multiple_align_functions(f, time, mu, lambda, option)
% Group-wise function alignment to specified mean
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
% option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS,expBayes)
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

fprintf('\n lambda = %5.1f \n', lambda);

% check dimension of time
a = size(time,1);
if a == 1
    time = time.';
end

[M, N] = size(f);

if option.smooth == 1
    f = smooth_data(f, option.sparam);
end

if option.showplot == 1
    figure(1); clf;
    plot(time, f, 'linewidth', 1);
    title('Original data', 'fontsize', 16);
    pause(0.1);
end

%% Compute the q-function of the plot
q = f_to_srvf(f,time);

%% Compute the q-function of the plot
mq = f_to_srvf(mu,time);

fn = zeros(M,N);
qn = zeros(M,N);
gam = zeros(N,size(q,1));
gam_dev = zeros(N,size(q,1));
mcmcopts.iter = option.MaxItr;
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
if option.parallel == 1
    parfor k = 1:N
        if (option.method=='expBayes')
            out_e = pairwise_align_bayes(mu, f(:,k), time, mcmcopts);
            gam(k,:) = out_e.gamma;
        else
            gam(k,:) = optimum_reparam(mq,q(:,k),time,lambda,option.method,option.w, ...
                mf(1,r), f(1,k,1));
        end
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
        fn(:,k) = interp1(t, f(:,k,1), (time(end)-time(1)).*gam(k,:) + time(1))';
        qn(:,k) = f_to_srvf(fn(:,k),time);
    end
else
    for k = 1:N
        if (option.method=='expBayes')
            out_e = pairwise_align_bayes(mu, f(:,k), time, mcmcopts);
            gam(k,:) = out_e.gamma;
        else
            gam(k,:) = optimum_reparam(mq,q(:,k),time,lambda,option.method,option.w, ...
                mf(1,r), f(1,k,1));
        end
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
        fn(:,k) = interp1(t, f(:,k,1), (time(end)-time(1)).*gam(k,:) + time(1))';
        qn(:,k) = f_to_srvf(fn(:,k),time);
    end
end

gamI = SqrtMeanInverse(gam);

%% Aligned data & stats
q0 = q;
mean_f0 = mean(f, 2);
std_f0 = std(f, 0, 2);
mean_fn = mean(fn, 2);
std_fn = std(fn, 0, 2);
mqn = mq;
fmean = mean(f(1,:))+cumtrapz(time,mqn.*abs(mqn));

fgam = zeros(M,N);
for ii = 1:N
    fgam(:,ii) = interp1(time, fmean, (time(end)-time(1)).*gam(ii,:) + time(1));
end
var_fgam = var(fgam,[],2);

stats.orig_var = trapz(time,std_f0.^2);
stats.amp_var = trapz(time,std_fn.^2);
stats.phase_var = trapz(time,var_fgam);

if option.showplot == 1
    figure(2); clf;
    plot((0:M-1)/(M-1), gam, 'linewidth', 1);
    axis square;
    title('Warping functions', 'fontsize', 16);
    
    figure(3); clf;
    plot(time, fn, 'LineWidth',1);
    title(['Warped data, \lambda = ' num2str(lambda)], 'fontsize', 16);
    
    figure(4); clf;
    plot(time, mean_f0, 'b-', 'linewidth', 1); hold on;
    plot(time, mean_f0+std_f0, 'r-', 'linewidth', 1);
    plot(time, mean_f0-std_f0, 'g-', 'linewidth', 1);
    title('Original data: Mean \pm STD', 'fontsize', 16);
    
    figure(5); clf;
    plot(time, mean_fn, 'b-', 'linewidth', 1); hold on;
    plot(time, mean_fn+std_fn, 'r-', 'linewidth', 1);
    plot(time, mean_fn-std_fn, 'g-', 'linewidth', 1);
    title(['Warped data, \lambda = ' num2str(lambda) ': Mean \pm STD'], 'fontsize', 16);
    
    figure(6); clf;
    plot(time, fmean, 'g','LineWidth',1);
    title(['f_{mean}, \lambda = ' num2str(lambda)], 'fontsize', 16);
end

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end

out.f0 = f;
out.time = time;
out.fn = fn;
out.qn = qn;
out.q0 = q0;
out.fmean = fmean;
out.mqn = mqn;
out.gam = gam;
out.psi = [];
out.stats = stats;
out.qun = [];
out.lambda = lambda;
out.method = option.method;
out.gamI = gamI;
out.rsamps = false;
