classdef elastic_lpcr_regression
    %elastic_pcr_regression A class to provide a SRVF logistic PCR regression
    % -------------------------------------------------------------------------
    % This class provides elastic logistic pcr regression for functional 
    % data using the SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_lpcr_regression(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: label vector 
    %   time: time vector of length M
    %
    %
    % elastic_lpcr_regression Properties:
    %   f - (M,N) % matrix defining N functions of M samples
    %   y - response vector of length N (-1/1)
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
    % elastic_lpcr_regression Methods:
    %   elastic_lpcr_regression - class constructor
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
        warp_data % fdawarp with alignment data
        alpha % intercept
        b % coefficient vector
        Loss % logistic loss
        PC % probability of classification
        pca % pca of aligned functional data
        y_labels % predicted labels
        
    end
    
    methods
        function obj = elastic_lpcr_regression(f, y, time)
            %elastic_regression Construct an instance of this class
            % Input:
            %   f: (M,N): matrix defining N functions of M samples
            %   y: response vector
            %   time: time vector of length M
            obj.time = time(:);
            obj.f = f;
            obj.y = y(:);
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
            b0 = zeros(no+1, 1);
            obj.b = fminunc(@(b) logit_optim(b,Phi,obj.y),b0,options);
            
            % Compute the loss
            obj.Loss = logit_loss(obj.b,Phi,obj.y);

            obj.alpha = obj.b(1);
            obj.b = obj.b(2:end);
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
                y_pred = phi(obj.alpha + a*obj.b);
                obj.y_labels = ones(size(y_pred));
                obj.y_labels(y_pred < 0.5) = -1;
                
                if (isempty(newdata.y))
                    obj.PC = NaN;
                else
                    obj.PC = sum(newdata.y(:) == obj.y_labels)./length(obj.y_labels);
                end
            else
                y_pred = phi(obj.alpha + obj.pca.coef*obj.b);
                obj.y_labels = ones(size(y_pred));
                obj.y_labels(y_pred < 0.5) = -1;
                obj.PC = sum(obj.y == obj.y_labels)./length(obj.y_labels);
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
