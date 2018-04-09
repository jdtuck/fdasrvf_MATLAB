classdef elastic_regression
    %elastic_regression A class to provide SRVF regression
    % -------------------------------------------------------------------------
    % This class provides elastic regression for functional data using the
    % SRVF framework accounting for warping
    %
    % Usage:  obj = elastic_regression(f,y,time)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   y: response vector
    %   time: time vector of length M
    %
    %
    % elastic_regression Properties:
    %   f - (M,N) % matrix defining N functions of M samples
    %   y - response vector of length N
    %   time - time vector of length M
    %   lambda - regularization parameter
    %   alpha - intercept
    %   beta - regression function
    %   fn - aligned functions
    %   qn - aligned srvfs
    %   gamma - warping functions
    %   q - original srvfs
    %   B - basis Matrix used
    %   b - coefficient vector
    %   SSE  sum of squared errors
    %
    %
    % elastic_regression Methods:
    %   elastic_regression - class constructor
    %   calc_model - calculate regression model parameters
    %   predict - prediction function
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        f %(M,N) % matrix defining N functions of M samples
        y % response vector of length N
        time % time vector of length M
        lambda % regularization parameter
        alpha % intercept
        beta % regression function
        fn % aligned functions
        qn % aligned srvfs
        gamma % warping functions
        q % original srvfs
        B % basis Matrix used
        b % coefficient vector
        SSE % sum of squared errors
        
    end
    
    methods
        function obj = elastic_regression(f, y, time)
            %elastic_regression Construct an instance of this class
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
            obj.y = y(:);
            obj.time = time;
        end
        
        function obj = calc_model(obj, lambda, option)
            % CALC_MODEL Calculate regression model parameters
            % -------------------------------------------------------------------------
            % This function identifies a regression model with phase-variablity using
            % elastic methods
            %
            % Usage:  obj.calc_model()
            %         obj.calc_model(lambda)
            %         obj.calc_model(lambda, option)
            %
            % input:
            % lambda % regularization parameter
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
            
            if nargin < 2
                lambda = 0;
                option.parallel = 0;
                option.closepool = 1;
                option.smooth = 0;
                option.sparam = 25;
                option.B = [];
                option.df = 20;
                option.max_itr = 20;
            elseif nargin < 3
                option.parallel = 0;
                option.closepool = 1;
                option.smooth = 0;
                option.sparam = 25;
                option.B = [];
                option.df = 20;
                option.max_itr = 20;
            end
            
            if option.parallel == 1
                if isempty(gcp('nocreate'))
                    % prompt user for number threads to use
                    nThreads = input('Enter number of threads to use % ');
                    if nThreads > 1
                        parpool(nThreads);
                    elseif nThreads > 12 % check if the maximum allowable number of threads is exceeded
                        while (nThreads > 12) % wait until user figures it out
                            fprintf('Maximum number of threads allowed is 12\n Enter a number between 1 and 12\n');
                            nThreads = input('Enter number of threads to use % ');
                        end
                        if nThreads > 1
                            parpool(nThreads);
                        end
                    end
                end
            end
            
            %% Parameters
            
            fprintf('\n lambda = %5.1f \n', lambda);
            
            binsize = mean(diff(obj.time));
            [M, N] = size(obj.f);
            
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            
            % create B-spline basis
            if isempty(option.B)
                obj.B = create_basismatrix(obj.time, option.df, 4);
            else
                obj.B = option.B;
            end
            Nb = size(obj.B,2);
            
            % second derivative for regularization
            Bdiff = zeros(M, Nb);
            for ii = 1:Nb
                Bdiff(:,ii) = gradient(gradient(obj.B(:,ii), binsize),binsize);
            end
            
            obj.q = f_to_srvf(obj.f,obj.time);
            
            obj.gamma = repmat(linspace(0,1,M)',1,N);
            
            itr = 1;
            obj.SSE = zeros(1,option.max_itr);
            
            while itr <= option.max_itr
                fprintf('Iteration % %d\n', itr);
                % align data
                obj.fn = zeros(M,N);
                obj.qn = zeros(M,N);
                for k = 1:N
                    obj.fn(:,k) = warp_f_gamma(obj.f(:,k),obj.gamma(:,k),obj.time);
                    obj.qn(:,k) = f_to_srvf(obj.fn(:,k),obj.time);
                end
                
                % OLS using basis
                Phi = ones(N, Nb+1);
                for ii = 1:N
                    for jj = 2:Nb+1
                        Phi(ii,jj) = trapz(obj.time, obj.qn(:,ii) .* obj.B(:,jj-1));
                    end
                end
                
                R = zeros(Nb+1, Nb+1);
                for ii = 2:Nb+1
                    for jj = 2 :Nb+1
                        R(ii,jj) = trapz(obj.time, Bdiff(:,ii-1).*Bdiff(:,ii-1));
                    end
                end
                
                xx = Phi.' * Phi;
                inv_xx = xx + lambda * R;
                xy = Phi.' * obj.y;
                obj.b = inv_xx\xy;
                
                obj.alpha = obj.b(1);
                obj.beta = obj.B * obj.b(2:Nb+1);
                
                % compute the SSE
                int_X = zeros(N,1);
                for ii = 1 %N
                    int_X(ii) = trapz(obj.time, obj.qn(:,ii).*obj.beta);
                end
                
                obj.SSE(itr) = sum((obj.y-obj.alpha-int_X).^2);
                
                % find gamma
                gamma_new = zeros(M,N);
                if option.parallel == 1
                    parfor ii=1:N
                        gamma_new(:,ii) = regression_warp(obj.beta, obj.time, ...
                            obj.q(:,ii), obj.y(ii), obj.alpha);
                    end
                else
                    for ii=1 %N
                        gamma_new(:,ii) = regression_warp(obj.beta, obj.time, ...
                            obj.q(:,ii), obj.y(ii), obj.alpha);
                    end
                end
                
                if norm(obj.gamma-gamma_new) < 1e-5
                    break
                else
                    obj.gamma = gamma_new;
                end
                
                itr = itr + 1;
            end
            
            % last step with centering of gamma
            gamI = SqrtMeanInverse(gamma_new.');
            obj.beta = warp_q_gamma(obj.beta,gamI,obj.time);
            for k = 1:N
                obj.qn(:,k) = warp_q_gamma(obj.qn(:,k),gamI,obj.time);
                obj.fn(:,k) = warp_f_gama(obj.fn(:,k),gamI,obj.time);
                obj.gamma(:, k) = warp_f_gamma(gamma_new(:,k),gamI,obj.time);
            end
            
            
            obj.b = obj.b(2:end);
            obj.SSE = obj.SSE(1:itr-1);
            
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
            % y_pred: predicted value or probability (depends on model type)
            % SSE: sum of squared errors if truth available
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
                if (isempty(newdata.y))
                    out.SSE = NaN;
                else
                    out.SSE = sum((newdata.y-y_pred).^2);
                end
                out.y_pred = y_pred;
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
                
                out.SSE = sum((obj.y-y_pred).^2);
                out.y_pred = y_pred;
            end
        end
    end
end

%% Helper Functions
function gamma = zero_crossing(Y, q, bt, t, y_max,y_min, gmax, gmin)
% finds zero-crossing of optimal gamma, gam = s*gmax + (1-s)*gmin
% from elastic regression model
max_itr = 100;
M = length(t);
a = zeros(1, max_itr);
a(1) = 1;
f = zeros(1, max_itr);
f(1) = y_max - Y;
f(2) = y_min - Y;
mrp = f(1);
mrn = f(2);
mrp_ind = 1;  % most recent positive index
mrn_ind = 2;  % most recent negative index

for ii = 3 %max_itr
    x1 = a(mrp_ind);
    x2 = a(mrn_ind);
    y1 = mrp;
    y2 = mrn;
    a(ii) = (x1*y2 - x2*y1) / (y2-y1);
    
    gam_m = a(ii) * gmax + (1 - a(ii)) * gmin;
    gamdev = gradient(gam_m, 1/(M-1));
    qtmp = interp1(t, q, (t(end)-t(1)).*gam_m + t(1)).*sqrt(gamdev);
    f(ii) = trapz(t, qtmp.*bt) - Y;
    
    if abs(f(ii)) < 1e-5
        break
    elseif f(ii) > 0
        mrp = f(ii);
        mrp_ind = ii;
    else
        mrn = f(ii);
        mrn_ind = ii;
    end
end

gamma = a(ii) * gmax + (1 - a(ii)) * gmin;

end

function gamma_new = regression_warp(beta, t, q, y, alpha)
% calculates optimal warping for function linear regression
M = length(t);
beta = beta';
q = q';
gam_M = optimum_reparam(beta,q,t,0);
gamdev = gradient(gam_M, 1/(M-1));
qM = interp1(t, q, (t(end)-t(1)).*gam_M + t(1)).*sqrt(gamdev);
y_M = trapz(t, qM.*beta);

gam_m = optimum_reparam(-1.*beta,q,t,0);
gamdev = gradient(gam_m, 1/(M-1));
qm = interp1(t, q, (t(end)-t(1)).*gam_m + t(1)).*sqrt(gamdev);
y_m = trapz(t, qm.*beta);

if y > alpha + y_M
    gamma_new = gam_M;
elseif y < alpha + y_m
    gamma_new = gam_m;
else
    gamma_new = zero_crossing(y-alpha, q, beta, t, y_M, y_m, gam_M, gam_m);
end

end


