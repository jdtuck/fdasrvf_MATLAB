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
        Loss % sum of squared errors
        pca % pca of aligned functional data
        n_classes % number of classes
        ylabels % predicted labels
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
            m = max(y);
            obj.n_classes = m;
            obj.Y = zeros(N,m);
            for ii =1:N
                obj.Y(ii,y(ii)) = 1;
            end
        end
        
        function obj = calc_model(obj, method, no, option)
            % CALC_MODEL Calculate regression model parameters
            % -------------------------------------------------------------------------
            % This function identifies a regression model with phase-variablity using
            % elastic methods
            %
            % Usage:  obj.calc_model(method, no, option)
            %         obj.calc_model(method, no)
            %
            % input:
            % method: string specifing pca method (options = "combined",
            %   "vert", or "horiz", default = "combined")
            % no: number of principal components
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
            method = lower(method);
            m = obj.n_classes;
            
            if nargin < 3
                option.parallel = 1;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 0;
                option.method = 'DP1';
                option.MaxItr = 20;
            end
            
            if option.smooth
                obj.f = smooth_data(obj.f,option.sparam);
                option.smooth = 0;
            end
            
            %% Align Data
            obj.warp_data = fdawarp(obj.f,obj.time);
            obj.warp_data = obj.warp_data.time_warping(0, option);
            
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
            out_pca = out_pca.calc_fpca(no);
            
            % LS using PCA basis
            Phi = ones(N1,no+1);
            Phi(:,2:(no+1)) = out_pca.coef;
            % find alpha and beta using bfgs
            options.Method = 'lbfgs';
            options.Display = 'off';
            b0 = zeros(m*(no+1), 1);
            obj.b = minFunc(@mlogit_optim,b0,options,Phi,obj.Y);
            
            % Compute the loss
            obj.LL = mlogit_loss(obj.b,Phi,obj.Y);
            
            B0 = reshape(obj.b, no+1, m);
            obj.alpha = B0(1,:);
            obj.b = B0(2:end,:);
            
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
            % PC: probability of classficiation if truth available
            omethod = obj.warp_data.method;
            lambda = obj.warp_data.lambda;
            m = obj.n_classes;
            M = length(obj.time);
            if (nargin>1)
                if (newdata.smooth)
                    newdata.f = smooth_data(newdata.f,newdata.sparam);
                end
                q1 = f_to_srvf(newdata.f,newdata.time);
                n = size(q1,2);
                y_pred = zeros(n,m);
                mq = obj.warp_data.mqn;
                fn = zeros(M,n);
                qn = zeros(M,n);
                gam = zeros(M,n);
                for ii = 1:n
                    gam(:,ii) = optimum_reparam(mq,q1(:,ii),obj.time,lambda,omethod);
                    fn(:,ii) = warp_f_gamma(newdata.f(:,ii),gam(:,ii),obj.time);
                    qn(:,ii) = f_to_srvf(fn(:,ii),obj.time);
                end
                m_new = sign(fn(obj.pca.id,:)).*sqrt(abs(fn(obj.pca.id,:)));
                qn1 = [qn; m_new];
                U = obj.pca.U;
                no = size(U,2);
                
                switch class(obj.pca)
                    case 'fdajpca'
                        C = obj.pca.C;
                        TT = length(obj.time);
                        mu_g = obj.pca.mu_g;
                        mu_psi = obj.pca.mu_psi;
                        vec = zeros(M,n);
                        psi = zeros(TT,n);
                        binsize = mean(diff(obj.time));
                        for i = 1:n
                            psi(:,i) = sqrt(gradient(gam(:,i),binsize));
                        end
                        
                        for i = 1:n
                            vec(:,i) = inv_exp_map(mu_psi, psi(:,i));
                        end
                        
                        g = [qn1; C*vec];
                        a = zeros(n,no);
                        for i = 1:n
                            for j = 1:no
                                a(i,j) = (g(:,i)-mu_g)*U(:,j);
                            end
                        end
                        
                    case 'fdavpca'
                        a = matrix(0,n,no);
                        for k = 1:no
                            for i = 1:n
                                a(i,k) = (qn1(:,i)-obj.pca.mqn)*U(:,k);
                            end
                        end
                        
                    case 'fdahpca'
                        a = zeros(n,no);
                        mu_psi = model.pca.mu;
                        vec = zeros(M,n);
                        TT = length(obj.time);
                        psi = zeros(TT,n);
                        binsize = mean(diff(objtime));
                        for i = 1:n
                            psi(:,i) = sqrt(gradient(gam(:,i),binsize));
                        end
                        
                        for i = 1:n
                            vec(:,i) = inv_exp_map(mu_psi, psi(:,i));
                        end
                        
                        vm = mean(obj.pca.vec,2);
                        
                        for k = 1:no
                            for i = 1:n
                                a(i,k) = sum((vec(:,i)-vm).*U(:,k));
                            end
                        end
                        
                    otherwise
                        error('invalid pca class');
                end
                
                for ii = 1:n
                    for jj = 1:m
                        y_pred(ii,jj) = obj.alpha(jj) + sum(a(ii,:).*obj.b(:,jj));
                    end
                end
                
                if (isempty(newdata.y))
                    y_pred = phi(reshape(y_pred,1,n*m));
                    y_pred = reshape(y_pred,n,m);
                    [~, obj.y_labels] = max(y_pred,[],2);
                    obj.PC = NaN;
                else
                    y_pred = phi(reshape(y_pred,1,n*m));
                    y_pred = reshape(y_pred,n,m);
                    [~, obj.y_labels] = max(y_pred,[],2);
                    obj.PC = zeros(1,m);
                    cls_set = 1:m;
                    for ii = 1:m
                        cls_sub = setdiff(cls_set,ii);
                        TP = sum(obj.y(obj.y_labels == ii) == ii);
                        FP = sum(obj.y(ismember(obj.y_labels,cls_sub)) == ii);
                        TN = sum(obj.y(ismember(obj.y_labels,cls_sub)) == ...
                            y_labels(ismember(obj.y_labels,cls_sub)));
                        FN = sum(ismember(obj.y(obj.y_labels==ii), cls_sub));
                        obj.PC(ii) = (TP+TN)/(TP+FP+FN+TN);
                    end
                    obj.PC = sum(option.y == y_labels)./length(y_labels);
                end
            else
                n = size(obj.pca.coef,1);
                y_pred = zeros(n,m);
                for ii = 1:n
                    for jj = 1:m
                        y_pred(ii,jj) = obj.alpha(jj) + sum(obj.pca.coef(ii,:).*obj.b(:,jj));
                    end
                end
                
                y_pred = phi(reshape(y_pred,1,n*m));
                y_pred = reshape(y_pred,n,m);
                [~, obj.y_labels] = max(y_pred,[],2);
                obj.PC = zeros(1,m);
                cls_set = 1:m;
                for ii = 1:m
                    cls_sub = setdiff(cls_set,ii);
                    TP = sum(obj.y(obj.y_labels == ii) == ii);
                    FP = sum(obj.y(ismember(obj.y_labels,cls_sub)) == ii);
                    TN = sum(obj.y(ismember(obj.y_labels,cls_sub)) == ...
                        y_labels(ismember(obj.y_labels,cls_sub)));
                    FN = sum(ismember(obj.y(obj.y_labels==ii), cls_sub));
                    obj.PC(ii) = (TP+TN)/(TP+FP+FN+TN);
                end
                obj.PC = sum(option.y == y_labels)./length(y_labels);
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

