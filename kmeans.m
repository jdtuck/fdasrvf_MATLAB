function out = kmeans(f,t,K,seeds,lambda,option)
% KMEANS K-Means clustering and alignment
% -------------------------------------------------------------------------
% This function clusters functions and aligns using the elastic square-root
% slope (srsf) framework.
%
% Usage:  out = kmeans(f,t,K)
%         out = kmeans(f,t,K,seeds)
%         out = kmeans(f,t,K,seeds,lambda)
%         out = kmeans(f,t,K,seeds,lambda,option)
%
% Input:
% f (M,N): matrix defining N functions of M samples
% t: time vector of length M
% K: number of clusters
% seeds: indexes of cluster center functions (default [])
% lambda: regularization parameter
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.alignment = true; % performa alignment
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting
% option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
% option.w = 0.0; % BFGS weight
% option.MaxItr = 20;  % maximum iterations
% option.thresh = 0.01; % cost function threshold
%
% Output:
% structure containing
% f0: original functions
% time: time
% fn: aligned functions in cell of clusters
% qn: aligned srvfs in cell of clusters
% gam: warping functions in cell of clusters
% q0: original srvfs
% labels: cluster labels
% templates: cluster centers
% templates_q: cluster srvf centers
% qun: cost function

if nargin < 4
    seeds = [];
    lambda = 0;
    option.parallel = 0;
    option.alignment = true;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
    option.thresh = 0.01;
elseif nargin < 5
    lambda = 0;
    option.parallel = 0;
    option.alignment = true;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
    option.thresh = 0.01;
elseif nargin < 6
    option.parallel = 0;
    option.alignment = true;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
    option.thresh = 0.01;
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

[M, N] = size(f);

if (isempty(seeds))
    template_ind = randsample(N,K);
else
    template_ind = seeds;
end
templates = zeros(M,K);
for i=1:K
    templates(:,i) = f(:,template_ind(i));
end
cluster_id = zeros(1,N);
qun = zeros(0, option.MaxItr);

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
templates_q = zeros(M,K);
for i=1:K
    templates_q(:,i) = q(:,template_ind(i));
end

%% Set initial using the original f space
fprintf('\nInitializing...\n');
f_temp = zeros(length(t),N);
q_temp = zeros(length(t),N);
for r = 1:option.MaxItr
    gam = {};
    Dy = zeros(K,N);
    qn = {};
    fn = {};
    fprintf('updating step: r=%d\n', r);
    if r == option.MaxItr
        fprintf('maximal number of iterations is reached. \n');
    end
    for i=1:K
        gamt = zeros(N,size(q,1));
        if option.parallel == 1
            parfor k = 1:N
                q_c = q(:,k); mq_c = templates_q(:,i);
                if (option.alignment)
                    gamt(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        templates(1,i), f(1,k));
                else
                    gamt(k,:) = linspace(0,1,M);
                end
                f_temp(:,k) = interp1(t, f(:,k), (t(end)-t(1)).*gamt(k,:) + t(1))';
                q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
                Dy(i,k) = sqrt(sum(trapz(t,(q_temp(:,k)-templates_q(:,i)).^2)));
            end
        else
            for k = 1:N
                q_c = q(:,k); mq_c = templates_q(:,i);
                if (option.alignment)
                    gamt(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        templates(1,i), f(1,k));
                else
                    gamt(k,:) = linspace(0,1,M);
                end
                f_temp(:,k) = interp1(t, f(:,k), (t(end)-t(1)).*gamt(k,:) + t(1))';
                q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
                Dy(i,k) = sqrt(sum(trapz(t,(q_temp(:,k)-templates_q(:,i)).^2)));
            end
        end
        gam{i} = gamt;
        qn{i} = q_temp;
        fn{i} = f_temp;
    end
    [~,cluster_id] = min(Dy,[],1);
    
    %% Normalization
    for i=1:K
        id = cluster_id == i;
        ftmp = fn{i}(:,id);
        gamtmp = gam{i}(id,:);
        gamI = SqrtMeanInverse(gamtmp);
        N1 = size(ftmp,2);
        f_temp1 = zeros(M,N1);
        q_temp1 = zeros(M,N1);
        gam1 = zeros(N1,M);
        if option.parallel == 1
            parfor k = 1:N1
                f_temp1(:,k) = interp1(t, ftmp(:,k), (t(end)-t(1)).*gamI + t(1))';
                q_temp1(:,k) = f_to_srvf(f_temp1(:,k),t);
                gam1(k,:) = interp1(t, gamtmp(k,:), (t(end)-t(1)).*gamI + t(1))';
            end
        else
            for k = 1:N1
                f_temp1(:,k) = interp1(t, ftmp(:,k), (t(end)-t(1)).*gamI + t(1))';
                q_temp1(:,k) = f_to_srvf(f_temp1(:,k),t);
                gam1(k,:) = interp1(t, gamtmp(k,:), (t(end)-t(1)).*gamI + t(1))';
            end
        end
        gam{i}(id,:) = gam1;
        qn{i}(:,id) = q_temp1;
        fn{i}(:,id) = f_temp1;
    end
    
    %% Template Identification
    qun_t = zeros(1,K);
    for i = 1:K
        id = cluster_id == i;
        old_templates_q = templates_q;
        templates_q(:,i) = mean(qn{i}(:,id),2);
        templates(:,i) = mean(fn{i}(:,id),2);
        qun_t(i) = norm(templates_q(:,i)-old_templates_q(:,i))/norm(old_templates_q(:,i));
    end
    qun(r) = mean(qun_t);
    
    if qun(r) < option.thresh || r >= option.MaxItr
        break;
    end
end

%% Output
ftmp = {};
qtmp = {};
gamtmp = {};
for i = 1:K
    id = cluster_id == i;
    ftmp{i} = fn{i}(:,id);
    qtmp{i} = qn{i}(:,id);
    gamtmp{i} = gam{i}(id,:);
end

if option.showplot == 1
    K = length(fn);
    colors = distinguishable_colors(K);
    figure(2); clf; hold all
    for k=1:K
        plot(t, templates(:,k), 'Color', colors(k,:));
    end
    title('Cluster Mean Functions');
    
    figure(3); clf;
    plot(t, ftmp{1}, 'Color', colors(1,:));
    hold all
    for k = 2:K
        plot(t, ftmp{k}, 'Color', colors(k,:));
    end
    for k = 1:K
        plot(t, templates(:,k), 'Color', colors(k,:));
    end
    title('Clustered Functions');
    
    figure(4); clf;
    plot(t, qtmp{1}, 'Color', colors(1,:));
    hold all
    for k = 2:K
        plot(t, qtmp{k}, 'Color', colors(k,:));
    end
    for k = 1:K
        plot(t, templates_q(:,k), 'Color', colors(k,:));
    end
    title('Clustered Functions');
end

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end

out.f0 = f;
out.time = t;
out.fn = ftmp;
out.qn = qtmp;
out.q0 = q;
out.labels = cluster_id;
out.templates = templates;
out.tempaltes_q = templates_q;
out.lambda = lambda;
out.method = option.method;
out.gamI = gamI;
out.qun = qun(1:r);
