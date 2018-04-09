classdef elastic_logistic
    %elastic_logistic A class to provide SRVF logistic regression
    % ---------------------------------------------------------------------
    % This class provides elastic logisitc regression for
    % functional data using the SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_logistic(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: response vector
    %   time: time vector of length M
    %
    %
    % elastic_logistic Properties:
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
    % elastic_logistic Methods:
    %   elastic_logistic - class constructor
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
        function obj = elastic_logistic(f, y, time)
            %elastic_logistic Construct an instance of this class
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
            % CALC_MODEL Calculate logistic regression model parameters
            % -------------------------------------------------------------------------
            % This function identifies a regression model with
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
            
            binsize = mean(diff(obj.time));
            [M, N] = size(obj.f);
            
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
            
            obj.q = f_to_srvf(obj.f,t);
            
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
                        Phi(ii,jj) = trapz(obj.time, obj.qn(:,ii) .* obj.B(:,jj-1));
                    end
                end
                
                % find alpha and beta using bfgs
                options.Method = 'lbfgs';
                options.Display = 'off';
                b0 = zeros(Nb+1, 1);
                obj.b = minFunc(@logit_optim,b0,options,Phi,y);
                
                obj.alpha = obj.b(1);
                obj.beta = obj.B * obj.b(2:Nb+1);
                
                % compute the lostic loss
                LL(itr) = logit_loss(obj.b, Phi, obj.y);
                
                % find gamma
                gamma_new = zeros(M,N);
                if option.parallel == 1
                    parfor ii=1:N
                        gamma_new(:,ii) = logistic_warp(obj.beta, obj.time, ...
                            obj.q(:,ii), obj.y(ii));
                    end
                else
                    for ii=1:N
                        gamma_new(:,ii) = logistic_warp(obj.beta, obj.time,  ...
                            obj.q(:,ii), obj.y(ii));
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
            if (exist(newdata))
                q1 = f_to_srvf(newdata.f,newdata.time);
                n = size(q1,2);
                y_pred = zeros(n,1);
                for ii = 1:n
                    difference = obj.q - repmat(q1(:,ii),1,size(obj.q,2));
                    dist = sum(abs(difference).^2).^(1/2);
                    [~, argmin] = min(dist);
                    q_tmp = warp_q_gamma(q1(:,ii), obj.gamma(:,argmin), newdata.time);
                    
                    y_pred(ii) = obj.alpha + trapz(newdata.time, q_tmp.' .* obj.beta);
                end
                y_pred = phi(y_pred);
                out.y_labels = ones(1,n);
                out.y_labels(y_pred < 0.5) = -1;
                if (isempty(newdata.y))
                    out.PC = NaN;
                else
                    TP = sum(newdata.y(out.y_labels == 1) == 1);
                    FP = sum(newdata.y(out.y_labels == -1) == 1);
                    TN = sum(newdata.y(out.y_labels == -1) == -1);
                    FN = sum(newdata.y(out.y_labels == 1) == -1);
                    out.PC = (TP+TN)/(TP+FP+FN+TN);
                end
                
            else
                n = size(obj.q,2);
                y_pred = zeros(n,1);
                for ii = 1:n
                    difference = obj.q - repmat(obj.q(:,ii),1,size(obj.q,2));
                    dist = sum(abs(difference).^2).^(1/2);
                    [~, argmin] = min(dist);
                    q_tmp = warp_q_gamma(obj.q(:,ii), obj.gamma(:,argmin), newdata.time);
                    
                    y_pred(ii) = obj.alpha + trapz(obj.time, q_tmp.' .* obj.beta);
                end
                y_pred = phi(y_pred);
                out.y_labels = ones(1,n);
                out.y_labels(y_pred < 0.5) = -1;
                TP = sum(obj.y(out.y_labels == 1) == 1);
                FP = sum(obj.y(out.y_labels == -1) == 1);
                TN = sum(obj.y(out.y_labels == -1) == -1);
                FN = sum(obj.y(out.y_labels == 1) == -1);
                out.PC = (TP+TN)/(TP+FP+FN+TN);
                
            end
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

