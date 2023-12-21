function out = align_fPCA(f,t,num_comp,lambda,option)
% ALIGN_FPCA Group-wise function alignment and PCA Extractions
% -------------------------------------------------------------------------
%
% This function aligns a collection of functions while extracting principal
% components.
%
% Usage:  out = align_fPCA(f,t)
%         out = align_fPCA(f,t,num_comp,lambda)
%         out = align_fPCA(f,t,num_comp,lambda,option)
%
% Input:
% f (M,N): matrix defining N functions of M samples
% t : time vector of length M
% num_comp: number of principal components to extract 9default = 3
% lambda: regularization parameter
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting
% option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
% option.MaxItr = 20;  % maximum iterations
%
% Output:
% structure containing
% fn: aligned functions
% qn: aligned srvfs
% q0: original srvfs
% fmean: function mean
% mqn: mean srvf
% gam: warping functions
% psi: srvf of gam
% vfpca: vertical fPCA structure
%   q_pca: srvf principal directions
%   f_pca: f principal directions
%   latent: latent values
%   coef: coefficients
%   U: eigenvectors
%   id: point used for f(0)
% stats: structure of statistics of alignment

if nargin < 3
    num_comp = 3;
    lambda = 0;
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.MaxItr = 20;
elseif nargin < 4
    lambda = 0;
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.MaxItr = 20;
elseif nargin < 5
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
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
coef = -2:2;
Nstd = length(coef);
NP = 1:num_comp;  % number of principal components
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

%% PCA step
fprintf('\nInitializing...\n');
mnq = mean(q,2);
dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
[~, min_ind] = min(dqq);
mq = q(:,min_ind);
mf = f(:,min_ind);
qhat_cent = q-mq*ones(1,N);
K = 1/M * (qhat_cent * qhat_cent.');
[U,~,~] = svd(K);

alpha_i = zeros(num_comp,N);
for ii = 1:num_comp
    for jj = 1:N
        alpha_i(ii,jj) = trapz(t,qhat_cent(:,jj).*U(:,ii));
    end
end

tmp = U(:,1:num_comp)*alpha_i;
if option.smooth==1
    mq = mq/norm(mq);
    for ii = 1:N
        if (sum(tmp(:,ii) ~=0))
            tmp(:,ii) = tmp(:,ii)/norm(tmp(:,ii));
        end
    end
end
qhat = mq*ones(1,N) + tmp;

gam = zeros(N,size(q,1));
if option.parallel == 1
    parfor k = 1:N
        gam(k,:) = optimum_reparam(qhat(:,k),q(:,k),t,lambda,option.method, ...
            mf(1), f(1,k));
    end
else
    for k = 1:N
        gam(k,:) = optimum_reparam(qhat(:,k),q(:,k),t,lambda,option.method, ...
            mf(1), f(1,k));
    end
end
%% Compute Mean
fprintf('Aligning %d functions in SRVF space to %d fPCA components...\n',N, num_comp);
MaxItr = option.MaxItr;
Dx = zeros(1,MaxItr);
f_temp = zeros(length(t),N);
q_temp = zeros(length(t),N);
for r = 1:MaxItr
    fprintf('updating step: r=%d\n', r);
    if r == MaxItr
        fprintf('maximal number of iterations is reached. \n');
    end
    
    % PCA Step
    for k=1:N
        f_temp(:,k) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
        q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
    end
    q(:,:,r+1) = q_temp;
    f(:,:,r+1) = f_temp;
    mq(:,r+1) = mean(q(:,:,r+1),2);
    mf(:,r+1) = mean(f(:,:,r+1),2);
    
    K = cov(q_temp.');
    [U,~,~] = svd(K);
    
    qhat_cent = q_temp - mq(:,r+1)*ones(1,N);
    alpha_i = zeros(num_comp,N);
    for ii = 1:num_comp
        for jj = 1:N
            alpha_i(ii,jj) = trapz(t,qhat_cent(:,jj).*U(:,ii));
        end
    end
    
    tmp = U(:,1:num_comp)*alpha_i;
    if option.smooth==1
        mq = mq/norm(mq);
        for ii = 1:N
            if (sum(tmp(:,ii) ~=0))
                tmp(:,ii) = tmp(:,ii)/norm(tmp(:,ii));
            end
        end
    end
    qhat = mq(:,r+1)*ones(1,N) + tmp;
    
    % Matching Step
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    gam = zeros(N,size(q,1));
    if option.parallel == 1
        parfor k = 1:N
            gam(k,:) = optimum_reparam(qhat(:,k),q(:,k,1),t,lambda,option.method, ...
                mf(1,r+1), f(1,k,r+1));
        end
    else
        for k = 1:N
            gam(k,:) = optimum_reparam(qhat(:,k),q(:,k,1),t,lambda,option.method, ...
                mf(1,r+1), f(1,k,r+1));
        end
    end
    
    psi = sqrt(diff(gam.').*(M-1));
    mu = ones(1,M-1);
    Dx1 = zeros(1,N);
    for ii=1:N
        Dx1(ii) = real(acos(sum(mu.'.*psi(:,ii)./(M-1))));
    end
    Dx(r+1) = max(Dx1);
    
    if abs(Dx(r+1)-Dx(r)) < 1e-4 || r >= MaxItr
        break;
    end
end

% last step with centering of gam
r = r+1;
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

% Compte final PCA in the q domain
mq_new = mean(qn,2);
id = round(length(t)/2);
m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));  % scaled version
mqn = [mq_new; mean(m_new)];
K = cov([qn;m_new]');

[U,S,~] = svd(K);
s = diag(S);
stdS = sqrt(s);

% compute the PCA in the q domain
q_pca = zeros(length(mq_new)+1,Nstd,num_comp);
for k = NP
    for i = 1:Nstd
        q_pca(:,i,k) = [mq_new; mean(m_new)] + coef(i)*stdS(k)*U(:,k);
    end
end

% compute the correspondence to the original function domain
f_pca = zeros(length(mq_new),Nstd,num_comp);
for k = NP
    for i = 1:Nstd
        f_pca(:,i,k) = cumtrapzmid(t,q_pca(1:end-1,i,k).*abs(q_pca(1:end-1,i,k)),sign(q_pca(end,i,k)).*(q_pca(end,i,k).^2), id);
    end
    fbar = mean(fn,2);
    fsbar = mean(f_pca(:,:,k),2);
    err = repmat(fbar-fsbar,1,Nstd);
    f_pca(:,:,k) = f_pca(:,:,k) + err;
end

% coefficients
c = zeros(size(qn,2),num_comp);
for jj = NP
    for ii = 1:size(qn,2)
        c(ii,jj) = dot([qn(:,ii);m_new(ii)]-mqn,U(:,jj));
    end
end

vfpca.q_pca = q_pca;
vfpca.f_pca = f_pca;
vfpca.latent = s;
vfpca.coef = c;
vfpca.id = id;
vfpca.mqn = mqn;
vfpca.time = t;
vfpca.c = c;
vfpca.U = U;

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
    
    cl = 'rbgmc';
    [~, ~, p1] = size(q_pca);
    num_plot = ceil(p1/3);
    for ii = 1:num_plot
        figure;
        for k1 = 1:3
            k = k1+(ii-1)*3;
            subplot(2,3,k1);
            for i = 1:length(coef)
                plot(t, q_pca(1:end-1,i,k), cl(i), 'linewidth', 2); hold on;
            end
            title(['q domain: PD ' num2str(k)], 'fontsize', 14);
            subplot(2,3,k1+3);
            for i = 1:length(coef)
                plot(t, f_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
            end
            title(['f domain: PD ' num2str(k)], 'fontsize', 14);
        end
    end
    cumm_coef = 100*cumsum(s)./sum(s);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
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
out.Dx = Dx(1:r);
out.lambda = lambda;
out.method = option.method;
out.gamI = gamI;
out.rsamps = false;
out.vfpca = vfpca;
