classdef elastic_mlogistic
    %elastic_mlogistic A class to provide wSRVF multionomial logistic
    % regression
    % -------------------------------------------------------------------------
    % This class provides elastic multinomial logisitc regression for
    % functional data using the SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_mlogistic(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: response vector
    %   time: time vector of length M
    %
    %
    % elastic_mlogistic Properties:
    %   f - (M,N) % matrix defining N functions of M samples
    %   y - response vector of length N
    %   time - time vector of length M
    %   alpha - intercept
    %   beta - regression function
    %   fn - aligned functions
    %   qn - aligned srvfs
    %   gamma - warping functions
    %   q - original srvfs
    %   B - basis Matrix used
    %   b - coefficient vector
    %   Loss - logistic loss
    %   n_classes - number of classes
    %
    %
    % elastic_mlogistic Methods:
    %   elastic_mlogistic - class constructor
    %   calc_model - calculate regression model parameters
    %   predict - prediction function
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        f % (M,N): matrix defining N functions of M samples
        y % response vector
        time % time vector with M samples
        alpha % intercept
        beta % regression function
        fn % aligned functions
        qn % aligned srvfs
        gamma % warping functions
        q % original srvfs
        B % basis Matrix used
        b % coefficient vector
        Loss % logistic loss
        n_classes % number of classes
    end
    
    methods
        function obj = elastic_mlogistic(f, y, time)
            %elastic_mlogistic Construct an instance of this class
            % Input:
            %   f: (M,N): matrix defining N functions of M samples
            %   y: response vector
            %   time: time vector of length M
            error('function not working properly');
            a = size(time,1);
            if (a ~=1)
                time = time';
            end
            
            obj.f = f;
            obj.y = y;
            obj.time = time;
        end
        
        function obj = calc_model(obj, option)
            % CALC_MODEL Calculate multinomial logistic regression model parameters
            % -------------------------------------------------------------------------
            % This function identifies a multinomial regression model with
            % phase-variablity using elastic methods
            %
            % Usage:  obj.calc_model()
            %         obj.calc_model(option)
            %
            % Input:
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
            % elastic_mlogistic object
            
            if nargin < 1
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
            
            [M, N] = size(obj.f);
            
            % code labs
            m = max(obj.y);
            Y = zeros(N,m);
            for ii =1:N
                Y(ii,obj.y(ii)) = 1;
            end
            
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            
            % create B-spline basis
            if isempty(option.B)
                obj.B = create_basismatrix(t, option.df, 4);
            else
                obj.B = option.B;
            end
            Nb = size(obj.B,2);
            
            obj.q = f_to_srvf(obj.f,obj.time);
            
            obj.gamma = repmat(linspace(0,1,M)',1,N);
            
            %% Main Loop
            itr = 1;
            LL = zeros(1,option.max_itr);
            while itr <= option.max_itr
                fprintf('Iteration: %d\n', itr);
                % align data
                obj.fn = zeros(M,N);
                obj.qn = zeros(M,N);
                for k = 1:N
                    obj.fn(:,k) = warp_f_gamma(obj.f(:,k),obj.gamma(:,k),obj.time);
                    obj.qn(:,k) = f_to_srvf(obj.fn(:,k),obj.time);
                end
                
                Phi = ones(N, Nb+1);
                for ii = 1:N
                    for jj = 2:Nb+1
                        Phi(ii,jj) = trapz(t, obj.qn(:,ii) .* obj.B(:,jj-1));
                    end
                end
                
                % find alpha and beta using bfgs
                options.Method = 'lbfgs';
                options.Display = 'off';
                b0 = zeros(m*(Nb+1), 1);
                obj.b = minFunc(@mlogit_optim,b0,options,Phi,Y);
                
                B0 = reshape(obj.b, Nb+1, m);
                obj.alpha = B0(1,:);
                obj.beta = zeros(M,m);
                for ii = 1:m
                    obj.beta(:,ii) = obj.B * B0(2:Nb+1,ii);
                end
                
                % compute the lostic loss
                LL(itr) = mlogit_loss(obj.b, Phi, Y);
                
                % find gamma
                gamma_new = zeros(M,N);
                if option.parallel == 1
                    parfor ii=1:N
                        gamma_new(:,ii) = mlogit_warp_grad(obj.alpha, obj.beta, ...
                            obj.time, obj.q(:,ii), Y(ii,:));
                    end
                else
                    for ii=1:N
                        gamma_new(:,ii) = mlogit_warp_grad(obj.alpha, obj.beta, ...
                            obj.time, obj.q(:,ii), Y(ii,:));
                    end
                end
                
                if norm(obj.gamma-gamma_new) < 1e-5
                    break
                else
                    obj.gamma = gamma_new;
                end
                
                itr = itr + 1;
            end
            obj.gamma = gamma_new;
            
            obj.b = obj.b(2:end);
            obj.Loss = LL(1:itr-1);
            obj.n_classes = m;
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
        end
        
        function out = predict(obj, newdata)
            % PREDICT Elastic Functional Regression Prediction
            % -------------------------------------------------------------------------
            % This function performs prediction on regression model on new
            % data if available or current stored data in object
            %
            % Usage:  obj.predict()
            %         obj.predict(newdata)
            %
            % Input:
            % newdata - struct containing new data for prediction
            % newdata.f - (M,N) matrix of functions
            % newdata.time - vector of time points
            % newdata.y - truth if available
            % newdata.smooth - smooth data if needed
            % newdata.sparam - number of times to run filter
            %
            % default options
            %
            % Output:
            % structure with fields:
            % y_labels: predicted labels
            % PC: probability of classification if truth available
            m = obj.n_classes;
            if (exist(newdata))
                q1 = f_to_srvf(newdata.f,newdata.time);
                n = size(q1,2);
                y_pred = zeros(n,m);
                
                for ii = 1:n
                    difference = obj.q - repmat(q1(:,ii),1,size(obj.q,2));
                    dist = sum(abs(difference).^2).^(1/2);
                    [~, argmin] = min(dist);
                    q_tmp = warp_q_gamma(q1(:,ii), obj.gamma(:,argmin), newdata.time);
                    
                    for jj = 1:m
                        y_pred(ii,jj) = obj.alpha(jj) + trapz(newdata.time, q_tmp.' .* obj.beta(:, jj));
                    end
                end
                y_pred = phi(reshape(y_pred,1,n*m));
                y_pred = reshape(y_pred,n,m);
                [~, out.y_labels] = max(y_pred,[],2);
                if (isempty(newdata.y))
                    out.PC = NaN;
                else
                    PC = zeros(1,m);
                    cls_set = 1:m;
                    for ii = 1:m
                        cls_sub = setdiff(cls_set,ii);
                        TP = sum(newdata.y(out.y_labels == ii) == ii);
                        FP = sum(newdata.y(ismember(out.y_labels,cls_sub)) == ii);
                        TN = sum(newdata.y(ismember(out.y_labels,cls_sub)) == ...
                            y_labels(ismember(out.y_labels,cls_sub)));
                        FN = sum(ismember(newdata.y(out.y_labels==ii), cls_sub));
                        PC(ii) = (TP+TN)/(TP+FP+FN+TN);
                    end
                    out.PC = sum(newdata.y == out.y_labels)./length(out.y_labels);
                end
            else
                n = size(obj.q,2);
                y_pred = zeros(n,m);
                for ii = 1:n
                    difference = obj.q - repmat(obj.q(:,ii),1,size(obj.q,2));
                    dist = sum(abs(difference).^2).^(1/2);
                    [~, argmin] = min(dist);
                    q_tmp = warp_q_gamma(obj.q(:,ii), obj.gamma(:,argmin), newdata.time);
                    
                    for jj = 1:m
                        y_pred(ii,jj) = obj.alpha(jj) + trapz(obj.time, q_tmp.' .* obj.beta(:, jj));
                    end
                end
                
                y_pred = phi(reshape(y_pred,1,n*m));
                y_pred = reshape(y_pred,n,m);
                [~, out.y_labels] = max(y_pred,[],2);
                PC = zeros(1,m);
                cls_set = 1:m;
                for ii = 1:m
                    cls_sub = setdiff(cls_set,ii);
                    TP = sum(obj.y(out.y_labels == ii) == ii);
                    FP = sum(obj.y(ismember(out.y_labels,cls_sub)) == ii);
                    TN = sum(obj.y(ismember(out.y_labels,cls_sub)) == ...
                        y_labels(ismember(out.y_labels,cls_sub)));
                    FN = sum(ismember(obj.y(out.y_labels==ii), cls_sub));
                    PC(ii) = (TP+TN)/(TP+FP+FN+TN);
                end
                out.PC = sum(obj.y == out.y_labels)./length(out.y_labels);
            end
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

