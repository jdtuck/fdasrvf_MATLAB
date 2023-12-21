classdef elastic_pcr_regression
    %elastic_pcr_regression A class to provide a SRVF PCR regression
    % -------------------------------------------------------------------------
    % This class provides elastic pcr regression for functional data using the
    % SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_pcr_regression(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: response vector
    %   time: time vector of length M
    %
    %
    % elastic_pcr_regression Properties:
    %   f - (M,N) % matrix defining N functions of M samples
    %   y - response vector of length N
    %   warp_data - fdawarp object of alignment
    %   pca - class dependent on fPCA method used object of fPCA
    %   information
    %   alpha - intercept
    %   b - coefficient vector
    %   SSE - sum of squared errors
    %
    %
    % elastic_pcr_regression Methods:
    %   elastic_pcr_regression - class constructor
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
        SSE % sum of squared errors
        pca % pca of aligned functional data
        y_pred % predicted values
    end
    
    methods
        function obj = elastic_pcr_regression(f, y, time)
            %elastic_regression Construct an instance of this class
            % Input:
            %   f: (M,N): matrix defining N functions of M samples
            %   y: response vector
            %   time: time vector of length M
            obj.time = time(:);
            obj.f = f;
            obj.y = y(:);
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
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
            % option.sparam = 25; % number of times to run filter
            % option.showplot = 0; % turns on and off plotting
            % option.max_itr = 20; % maximum number of iterations
            %
            % output %
            % elastic_regression object
            method = lower(method);
            
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
            
            % OLS using PCA basis
            lam = 0;
            R = 0;
            N = size(obj.f,2);
            Phi = ones(N,no+1);
            Phi(:,2:(no+1)) = out_pca.coef;
            xx = Phi.' * Phi;
            xy = Phi.' * obj.y;
            obj.b = (xx + lam * R)\xy;
            obj.alpha = obj.b(1);
            obj.b = obj.b(2:end);
            
            % compute the SSE
            int_X = zeros(N,1);
            for ii = 1:N
                int_X(ii) = out_pca.coef(ii,:) * obj.b;
            end
            
            obj.SSE = sum((obj.y-obj.alpha-int_X).^2);
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
            % y_pred: predicted value or probability (depends on model type)
            % SSE: sum of squared errors if truth available
            omethod = obj.warp_data.method;
            lambda = obj.warp_data.lambda;
            M = length(obj.time);
            if (nargin>1)
                if (newdata.smooth)
                    newdata.f = smooth_data(newdata.f,newdata.sparam);
                end
                q1 = f_to_srvf(newdata.f,newdata.time);
                n = size(q1,2);
                obj.y_pred = zeros(n,1);
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
                        mu_psi = obj.pca.mu;
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
                    obj.y_pred(ii) = obj.alpha + sum(a(ii,:).*obj.b);
                end
                
                if (isempty(newdata.y))
                    obj.SSE = NaN;
                else
                    obj.SSE = sum((newdata.y-obj.y_pred).^2);
                end
                obj.y_pred = obj.y_pred;
            else
                n = size(obj.pca.coef,1);
                obj.y_pred = zeros(n,1);
                for ii = 1:n
                    obj.y_pred(ii) = obj.alpha + obj.pca.coef(ii,:)*obj.b;
                end
                
                obj.SSE = sum((obj.y-obj.y_pred).^2);
                obj.y_pred = obj.y_pred;
            end
        end
    end
end