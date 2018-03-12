function out = elastic_mlogistic(f, y, t, option)
% This function identifies a multinomial regression model with 
% phase-variablity using elastic methods
%
% input:
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
% Loss: logistic loss
% n_classes: number of classes
% type: model type

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

% code labs 
m = max(y);
Y = zeros(N,m);
for ii =1:N
    Y(ii,y(ii)) = 1;
end

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
    b0 = zeros(m*(Nb+1), 1);
    b = minFunc(@mlogit_optim,b0,options,Phi,Y);
    
    B0 = reshape(b, Nb+1, m);
    alpha = B0(1,:);
    beta = zeros(M,m);
    for ii = 1:m
        beta(:,ii) = B * B0(2:Nb+1,ii);
    end
    
    % compute the lostic loss
    LL(itr) = mlogit_loss(b, Phi, Y);
        
    % find gamma
    gamma_new = zeros(M,N);
    if option.parallel == 1
        parfor ii=1:N
            gamma_new(:,ii) = mlogit_warp_grad(alpha, beta, t, q(:,ii), Y(ii,:));
        end
    else
        for ii=1:N
            gamma_new(:,ii) = mlogit_warp_grad(alpha, beta, t, q(:,ii), Y(ii,:));
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
out.Loss = LL(1:itr-1);
out.n_classes = m;
out.type = 'mlogistic';

if option.parallel == 1 && option.closepool == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end

end

%% Helper Functions

function nll = mlogit_loss(b, X, Y)
% calculates multinomial logistic loss (negative log-likelihood)
[N, m] = size(Y);
M = size(X,2);
B = reshape(b,M,m);
Yhat = X * B;
Yhat = Yhat - repmat(min(Yhat,[],2),1,m);
Yhat = exp(-1.*Yhat);
% l1-normalize
Yhat = Yhat./repmat(sum(Yhat,2),1,m);

Yhat = Yhat .* Y;
nll = sum(log(sum(Yhat,2)));
nll = nll ./ (-1*N);
end

function grad = mlogit_gradient(b, X, Y)
% calculates gradient of the multinomial logistic loss
[N, m] = size(Y);
M = size(X,2);
B = reshape(b,M,m);
Yhat = X * B;
Yhat = Yhat - repmat(min(Yhat,[],2),1,m);
Yhat = exp(-1.*Yhat);
% l1-normalize
Yhat = Yhat./repmat(sum(Yhat,2),1,m);

Yhat1 = Yhat .* Y;
Yhat1 = Yhat1./repmat(sum(Yhat1,2),1,m);
Yhat = Yhat - Yhat1;
grad = X.' * Yhat;
grad = grad/(-1*N);
grad = reshape(grad, M*m, 1);
end

function [nll, g] = mlogit_optim(b, X, Y)
% function for call to optimizer
nll = mlogit_loss(b, X, Y);
if nargout > 1
    g = mlogit_gradient(b, X, Y);
end
end

function gamma = mlogit_warp_grad(alpha, beta, t, q, y)
% calculates optimal warping for functional multinomial logistic regression
max_itr=int32(8000);
tol=1e-10;
delta=0.008;
display=int32(0);
m1 = length(t);
m = size(beta,2);
y = int32(y);

alpha = alpha./norm(alpha);
q = q./norm(q);
for ii = 1:m
    beta(:,ii) = beta(:,ii)./norm(beta(:,ii));
end
beta1 = zeros(1,m1*m);
for ii = 1:m
    beta1(((ii-1)*m1+1):(ii*m1)) = beta(:, ii);
end
gami = linspace(0,1,m1);

gamma = mlogit_warp(alpha, beta1, t, q, y, gami, max_itr, tol, delta, display); 

end

