function out_warp = gauss_model(out_warp, n, sort_samples)
% GAUSS_MODEL Gaussian gnerative model
% -------------------------------------------------------------------------
% This function models the functional data using a Gaussian model extracted
% from the principal components of the srvfs
%
% Usage: samples = gauss_model(out_warp, n, sort_samples)
%
% Inputs:
% fn: matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
% time: vector of size \eqn{N} describing the sample points
% qn: matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned srvfs
% gam: warping functions
% n: number of random samples (n = 1)
% sort_samples: sort samples (e.g., = F)
%
% Output:
% Structure containing
% fs: random aligned samples
% gams: random warping function samples
% ft: random function samples

fn = out_warp.fn;
time = out_warp.time;
qn = out_warp.qn;
gam = out_warp.gam;

if nargin < 5
    n = 1;
    sort_samples = false;
elseif nargin <6
    sort_samples = false;
end
%% Separated and Warped Data
% sampling from the estimated model
[M, ~] = size(fn);
binsize = mean(diff(time));

% compute mean and covariance in q-domain
id = round(length(time)/2);
q_new = qn;
mq_new = mean(qn,2);
m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));
mqn2 = [mq_new; mean(m_new)];
C = cov([q_new;m_new]');

q_s = mvnrnd(mqn2', C, n)';

% compute the correspondence to the original function domain
f_s = zeros(M,n);
for k = 1:n
    f_s(:,k) = cumtrapzmid(time,q_s(1:(end-1),k).*abs(q_s(1:(end-1),k)),sign(q_s(end,k))*(q_s(end,k)^2),id);
end
fbar = mean(fn,2);
fsbar = mean(f_s,2);
err = repmat(fbar-fsbar,1,n);
f_s = f_s + err;

% random warping generation
rgam = randomGamma(gam.',n);

% combine samples
if sort_samples
    %%%% sort functions and warpings
    mx = max(f_s);
    [~, seq1] = sort(mx);
    
    % compute the psi-function
    psi = zeros(n,size(rgam,2));
    len = zeros(1,n);
    ip = zeors(1,n);
    for i = 1:n
        psi(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
        ip(i) = ones(1,M)*psi(i,:)'/M;
        len(i) = acos(ones(1,M)*psi(i,:)'/M);
    end
    [~, seq2] = sort(len);
    
    % combine x-variability and y-variability
    f_c = zeros(size(fn,1),n);
    for k = 1:n
        f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,seq1(k)), invertGamma(rgam(seq2(k),:)')');
        while sum(isnan(f_c(:,k))) >= 1
            rgam2 = randomGamma(gam.',1);
            f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,k), invertGamma(rgam2'));
        end
    end
else
    % compute the psi-function
    psi = zeros(n,size(rgam,2));
    len = zeros(1,n);
    ip = zeros(1,n);
    for i = 1:n
        psi(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
        ip(i) = ones(1,M)*psi(i,:)'/M;
        len(i) = acos(ones(1,M)*psi(i,:)'/M);
    end
    
    % combine x-variability and y-variability
    f_c = zeros(size(fn,1),n);
    for k = 1:n
        f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,k), invertGamma(rgam(k,:)')');
        while sum(isnan(f_c(:,k))) >= 1
            rgam2 = randomGamma(gam.',1);
            f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,k), invertGamma(rgam2'));
        end
    end
end

out_warp.fs = f_s;
out_warp.gams = rgam;
out_warp.ft = rgam;
out_warp.qs = q_s(1:M,:);
out_warp.rsamps = true;
end

function rgam = randomGamma(gam,num)

[mu,~,vec] = SqrtMean(gam);

K = cov(vec);
[U,S,~] = svd(K);
Sig = diag(S);
n = 5;
T = length(vec);
time = linspace(0,1,T);
time = time(:);

vm = mean(vec);
rgam = zeros(num, length(gam));
for k=1:num
    
    a = randn(1,n);
    v = zeros(size(vm));
    for i=1:n
        v = v + a(i)*sqrt(Sig(i))*U(:,i)';
    end
    psi = exp_map(mu, v);
    
    gam0 = cumtrapz(time,psi.^2);
    rgam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
    
end
end