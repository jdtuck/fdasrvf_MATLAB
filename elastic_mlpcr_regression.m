classdef elastic_mlpcr_regression
    %elastic_mlpcr_regression A class to provide a SRVF multinomial logistic PCR regression
    % -------------------------------------------------------------------------
    % This class provides elastic multinomial logistic pcr regression for functional
    % data using the SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_mlpcr_regression(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: label vector
    %   time: time vector of length M
    %
    %
    % elastic_mlpcr_regression Properties:
    %   f - (M,N) % matrix defining N functions of M samples
    %   y - response vector of length N
    %   Y - coded label matrix
    %   warp_data - fdawarp object of alignment
    %   pca - class dependent on fPCA method used object of fPCA
    %   information
    %   alpha - intercept
    %   b - coefficient vector
    %   Loss - logistic loss
    %   PC - probability of classification
    %   ylabels - predicted labels
    %
    %
    % elastic_mlpcr_regression Methods:
    %   elastic_mlpcr_regression - class constructor
    %   calc_model - calculate regression model parameters
    %   predict - prediction function
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  18-Mar-2018
    
    properties
        f % (M,N): matrix defining N functions of M samples
        time % time vector of length M
        y % response vector of length N
        Y % coded label matrix
        warp_data % fdawarp with alignment data
        alpha % intercept
        b % coefficient vector
        Loss % multinomial logistic loss
        pca % pca of aligned functional data
        n_classes % number of classes
        y_labels % predicted labels
        PC % probability of classification
        
    end
    
    methods
        function obj = elastic_mlpcr_regression(f, y, time)
            %elastic_regression Construct an instance of this class
            % Input:
            %   f: (M,N): matrix defining N functions of M samples
            %   y: response vector
            %   time: time vector of length M
            obj.time = time(:);
            obj.f = f;
            obj.y = y(:);
            
            % code labels
            N = length(obj.y);
            m = max(y);
            obj.n_classes = m;
            obj.Y = zeros(N,m);
            for ii =1:N
                obj.Y(ii,y(ii)) = 1;
            end
        end
        
        function obj = calc_model(obj, method, var_exp, option)
            % CALC_MODEL Calculate regression model parameters
            % -------------------------------------------------------------------------
            % This function identifies a regression model with phase-variablity using
            % elastic methods
            %
            % Usage:  obj.calc_model(method, var_exp, option)
            %         obj.calc_model(method, var_exp)
            %
            % input:
            % method: string specifing pca method (options = "combined",
            %   "vert", or "horiz", default = "combined")
            % var_exp: compute number of pcs based on value percent variance explained (default = 0.99)
            % option: option for alignment
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
            % output %
            % elastic_regression object
            
            arguments
                obj
                method
                var_exp = 0.99;
                option.parallel = 1;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 0;
                option.method = 'DP1';
                option.MaxItr = 20;
            end

            method = lower(method);
            m = obj.n_classes;
            
            
            if option.smooth
                obj.f = smooth_data(obj.f,option.sparam);
                option.smooth = 0;
            end
            
            %% Align Data
            obj.warp_data = fdawarp(obj.f,obj.time);
            obj.warp_data = obj.warp_data.time_warping(0, ...
                parallel=option.parallel, closepool=option.closepool, ...
                smooth=option.smooth, sparam=option.sparam, ...
                method=option.method, MaxItr=option.MaxItr);
            
            switch method
                case 'combined'
                    out_pca = fdajpca(obj.warp_data);
                case 'vert'
                    out_pca = fdavpca(obj.warp_data);
                case 'horiz'
                    out_pca = fdahpca(obj.warp_data);
                otherwise
                    error('Invalid Method')
            end
            out_pca = out_pca.calc_fpca(var_exp);
            no = size(out_pca.coef,2);
            
            % LS using PCA basis
            Phi = ones(size(out_pca.coef,1),no+1);
            Phi(:,2:(no+1)) = out_pca.coef;
            % find alpha and beta using bfgs
            options = optimoptions("fminunc",Algorithm="quasi-newton", ...
                SpecifyObjectiveGradient=true,Display="off");
            b0 = zeros(m*(no+1), 1);
            obj.b = fminunc(@(b) mlogit_optim(b,Phi,obj.Y),b0,options);
            
            % Compute the loss
            obj.Loss = mlogit_loss(obj.b,Phi,obj.Y);
            
            B0 = reshape(obj.b, no+1, m);
            obj.alpha = B0(1,:);
            obj.b = B0(2:end,:);
            obj.pca = out_pca;
        end
        
        function obj = predict(obj, newdata)
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
            % newdata.time - vector of time points (must match the training grid)
            % newdata.y - truth if available
            % newdata.smooth - smooth data if needed
            % newdata.sparam - number of times to run filter
            %
            % default options
            %
            % Output:
            % structure with fields:
            % y_labels: predicted labels
            % PC: probability of classficiation if truth available
            if (nargin>1)
                if (newdata.smooth)
                    newdata.f = smooth_data(newdata.f,newdata.sparam);
                end
                % align and project the new data onto the fPCA basis the
                % model was fit on
                new_pca = obj.pca.project(newdata.f);
                a = new_pca.new_coef(:,1:size(obj.pca.coef,2));
                y_pred = softmax_rows(obj.alpha + a*obj.b);
                [~, obj.y_labels] = max(y_pred,[],2);
                
                if (isempty(newdata.y))
                    obj.PC = NaN;
                else
                    obj.PC = sum(newdata.y(:) == obj.y_labels)./length(obj.y_labels);
                end
            else
                y_pred = softmax_rows(obj.alpha + obj.pca.coef*obj.b);
                [~, obj.y_labels] = max(y_pred,[],2);
                obj.PC = sum(obj.y == obj.y_labels)./length(obj.y_labels);
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
% softmax, P(class j) proportional to exp(Yhat(:,j))
Yhat = softmax_rows(Yhat);

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
% softmax, P(class j) proportional to exp(Yhat(:,j))
Yhat = softmax_rows(Yhat);

grad = X.' * (Yhat - Y);
grad = grad/N;
grad = reshape(grad, M*m, 1);
end

function [nll, g] = mlogit_optim(b, X, Y)
% function for call to optimizer
nll = mlogit_loss(b, X, Y);
if nargout > 1
    g = mlogit_gradient(b, X, Y);
end
end

function P = softmax_rows(Z)
% row-wise softmax, P(i,j) = exp(Z(i,j)) / sum_k exp(Z(i,k))
P = exp(Z - max(Z,[],2));
P = P./sum(P,2);
end

