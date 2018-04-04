classdef fdawarp
    %fdawarp A class to provide a SRVF functional data analysis
    % -------------------------------------------------------------------------
    % This class provides alignment methods for functional data using the
    % SRVF framework
    %
    % Usage:  obj = fdawarp(f,t)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   time: time vector of length M
    %
    %
    % fdawarp Properties:
    %    f - (M,N): matrix defining N functions of M samples
    %    time - time vector of length M
    %    fn - aligned functions
    %    qn - aligned srvfs
    %    q0 - initial srvfs
    %    fmean - function mean
    %    mqn - mean srvf
    %    gam - warping functions
    %    psi - srvf of warping functions
    %    stats - alignment statistics
    %    qun - cost function
    %    lambda - lambda
    %    method - optimization method
    %    gamI - invserse warping function
    %    rsamps - random samples
    %    fs - random aligned functions
    %    gams - random warping functions
    %    ft - random warped functions
    %    qs - random aligned srvfs
    %    type - alignment type
    %    mcmc - mcmc output if bayesian
    %
    %
    % fdawarp Methods:
    %   fdawarp - class constructor
    %   time_warping - align functions and find karcher mean
    %   time_warping_median - align functions and find karcher median
    %   multiple_align_functions - align functions to a specified mean
    %   gauss_model - model functions using a generative gaussian model on
    %       fPCA space
    %   plot - plot results and functions in object
    %          
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        f      % (M,N): matrix defining N functions of M samples
        time   % time vector of length M
        fn     % aligned functions
        qn     % aligned srvfs
        q0     % initial srvfs
        fmean  % function mean
        mqn    % mean srvf
        gam    % warping functions
        psi    % srvf of warping functions
        stats  % alignment statistics
        qun    % cost function
        lambda % lambda
        method % optimization method
        gamI   % invserse warping function
        rsamps % random samples
        fs     % random aligned functions
        gams   % random warping functions
        ft     % random warped functions
        qs     % random aligned srvfs
        type   % alignment type
        mcmc   % mcmc output if bayesian
        
    end
    
    methods
        function obj = fdawarp(f,time)
            %fdawarp Construct an instance of this class
            % Input:
            %   f: (M,N): matrix defining N functions of M samples
            %   time: time vector of length M
            
            % check dimension of time
            a = size(time,1);
            if a == 1
                time = time';
            end
            
            obj.f = f;
            obj.time = time;
        end
        
        function obj = time_warping(obj,lambda,option)
            % TIME_WARPING Group-wise function alignment
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the elastic square-root
            % slope (srsf) framework.
            %
            % Usage:  obj.time_warping()
            %         obj.time_warping(lambda)
            %         obj.time_warping(lambda,option)
            %
            % Input:
            % lambda: regularization parameter
            %
            % default options
            % option.parallel = 0; % turns on MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
            % option.w = 0.0; % BFGS weight
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdawarp object
            if nargin < 2
                lambda = 0;
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
                option.MaxItr = 20;
            elseif nargin < 3
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
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
            
            binsize = mean(diff(obj.time));
            [M, N] = size(obj.f);
            
            f1 = obj.f;
            if option.smooth == 1
                f1 = smooth_data(f1, option.sparam);
            end
            
            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time);
            
            %% Set initial using the original f space
            fprintf('\nInitializing...\n');
            mnq = mean(q,2);
            dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
            [~, min_ind] = min(dqq);
            mq = q(:,min_ind);
            mf = obj.f(:,min_ind);
            
            gam_o = zeros(N,size(q,1));
            if option.parallel == 1
                parfor k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                        mf(1), f1(1,k,1));
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                        mf(1), f1(1,k,1));
                end
            end
            gamI_o = SqrtMeanInverse(gam_o);
            mf = warp_f_gamma(mf,gamI_o,obj.time);
            mq = f_to_srvf(mf,obj.time);
            
            %% Compute Mean
            fprintf('Computing Karcher mean of %d functions in SRVF space...\n',N);
            ds = inf;
            MaxItr = option.MaxItr;
            qun_o = zeros(1,MaxItr);
            f_temp = zeros(length(obj.time),N);
            q_temp = zeros(length(obj.time),N);
            for r = 1:MaxItr
                fprintf('updating step: r=%d\n', r);
                if r == MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end
                
                % Matching Step
                clear gam gam_dev;
                % use DP to find the optimal warping for each function w.r.t. the mean
                gam_o = zeros(N,size(q,1));
                gam_dev = zeros(N,size(q,1));
                if option.parallel == 1
                    parfor k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                            mf(1,r), f1(1,k,1));
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                    end
                else
                    for k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                            mf(1,r), f1(1,k,1));
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                    end
                end
                q(:,:,r+1) = q_temp;
                f1(:,:,r+1) = f_temp;
                
                ds(r+1) = sum(simps(obj.time, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
                    lambda*sum(simps(obj.time, (1-sqrt(gam_dev')).^2));
                
                % Minimization Step
                % compute the mean of the matched function
                mq(:,r+1) = mean(q(:,:,r+1),2);
                mf(:,r+1) = mean(f1(:,:,r+1),2);
                
                qun_o(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
                if qun_o(r) < 1e-2 || r >= MaxItr
                    break;
                end
            end
            
            % last step with centering of gam
            r = r+1;
            if option.parallel == 1
                parfor k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                        mf(1,r), f1(1,k,1));
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                        mf(1,r), f1(1,k,1));
                end
            end
            gamI_o = SqrtMeanInverse(gam_o);
            mq(:,r+1) = warp_q_gamma(mq(:,r),gamI_o,obj.time);
            for k = 1:N
                q(:,k,r+1) = warp_q_gamma(q(:,k,r),gamI_o,obj.time);
                f1(:,k,r+1) = warp_f_gamma(f1(:,k,r),gamI_o,obj.time);
                gam_o(k,:) = interp1(obj.time, gam_o(k,:), (obj.time(end)-obj.time(1)).*gamI_o + obj.time(1));
            end
            
            %% Aligned data & stats
            obj.fn = f1(:,:,r+1);
            obj.qn = q(:,:,r+1);
            obj.q0 = q(:,:,1);
            std_f0 = std(f1, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            obj.mqn = mq(:,r+1);
            obj.fmean = mean(obj.f(1,:))+cumtrapz(obj.time,obj.mqn.*abs(obj.mqn));
            
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean,gam_o(ii,:),obj.time);
            end
            var_fgam = var(fgam,[],2);
            
            obj.stats.orig_var = trapz(obj.time,std_f0.^2);
            obj.stats.amp_var = trapz(obj.time,std_fn.^2);
            obj.stats.phase_var = trapz(obj.time,var_fgam);
            
            obj.gam = gam_o.';
            [~,fy] = gradient(obj.gam,binsize,binsize);
            obj.psi = sqrt(fy+eps);
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            obj.qun = qun_o(1:r-1);
            obj.lambda = lambda;
            obj.method = option.method;
            obj.gamI = gamI_o;
            obj.rsamps = false;
            obj.type = 'mean';
        end
        
        function obj = time_warping_median(obj,lambda,option)
            % TIME_WARPING_MEDIAN Group-wise function alignment to median
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the elastic square-root
            % slope (srsf) framework to the median
            %
            % Usage:  obj.time_warping_median()
            %         obj.time_warping_median(lambda)
            %         obj.time_warping_median(lambda,option)
            %
            % Input:
            % lambda: regularization parameter
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, DP2, SIMUL, RBFGS)
            % option.w = 0.0; % BFGS weight
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdawarp object
            if nargin < 2
                lambda = 0;
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
                option.MaxItr = 20;
            elseif nargin < 3
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
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
            
            t = obj.time;
            
            binsize = mean(diff(t));
            [M, N] = size(obj.f);
            f1 = obj.f;
            
            if option.smooth == 1
                f1 = smooth_data(f1, option.sparam);
            end
            
            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,t);
            
            %% Set initial using the original f space
            fprintf('\nInitializing...\n');
            mnq = mean(q,2);
            dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
            [~, min_ind] = min(dqq);
            mq = q(:,min_ind);
            mf = f1(:,min_ind);
            
            gam_o = zeros(N,size(q,1));
            if option.parallel == 1
                parfor k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        mf(1), f1(1,k,1));
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        mf(1), f1(1,k,1));
                end
            end
            
            gamI_o = SqrtMeanInverse(gam_o);
            mf = warp_f_gamma(mf,gamI_o,t);
            mq = f_to_srvf(mf,t);
            
            %% Compute Mean
            fprintf('Computing Karcher median of %d functions in SRVF space...\n',N);
            ds = inf;
            MaxItr = option.MaxItr;
            qun_o = zeros(1,MaxItr);
            f_temp = zeros(length(t),N);
            q_temp = zeros(length(t),N);
            for r = 1:MaxItr
                fprintf('updating step: r=%d\n', r);
                if r == MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end
                
                % Matching Step
                clear gam gam_dev;
                % use DP to find the optimal warping for each function w.r.t. the mean
                gam_o = zeros(N,size(q,1));
                gam_dev = zeros(N,size(q,1));
                vtil = zeros(M,N);
                dtil = zeros(1,N);
                if option.parallel == 1
                    parfor k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                            mf(1,r), f1(1,k,1));
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),t);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
                        v = q_temp(:,k) - mq_c
                        d = sqrt(trapz(t, v.*v));
                        vtil(:,k) = v/d;
                        dtil(k) = 1/d;
                    end
                else
                    for k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                            mf(1,r), f1(1,k,1));
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),t);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t);
                        v = q_temp(:,k) - mq_c;
                        d = sqrt(trapz(t, v.*v));
                        vtil(:,k) = v/d;
                        dtil(k) = 1/d;
                    end
                end
                q(:,:,r+1) = q_temp;
                f1(:,:,r+1) = f_temp;
                
                ds_tmp = sqrt(sum(simps(t,(q(:,:,r+1)-mq(:,r)*ones(1,N)).^2)))+lambda*sum(simps(t,(1-sqrt(gam_dev')).^2));
                if (isreal(ds_tmp))
                    ds(r+1) = ds_tmp;
                else
                    ds(r+1) = abs(ds_tmp);
                end
                
                
                % Minimization Step
                % compute the mean of the matched function
                stp = .3;
                vbar = sum(vtil,2)*sum(dtil)^(-1);
                mq(:,r+1) = mq(:,r) + stp*vbar;
                mf(:,r+1) = median(f1(1,:,1)) + cumtrapz(t, mq(:,r+1).*abs(mq(:,r+1)));
                
                qun_o(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
                if qun_o(r) < 1e-2 || r >= MaxItr
                    break;
                end
            end
            
            % last step with centering of gam
            r = r+1;
            if option.parallel == 1
                parfor k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        mf(1,r), f1(1,k,1));
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method,option.w, ...
                        mf(1,r), f1(1,k,1));
                end
            end
            
            gamI_o = SqrtMeanInverse(gam_o);
            mq(:,r+1) = warp_q_gamma(mq(:,r),gamI_o,obj.time);
            for k = 1:N
                q(:,k,r+1) = warp_q_gamma(q(:,k,r),gamI_o,obj.time);
                f1(:,k,r+1) = warp_f_gamma(f1(:,k,r),gamI_o,obj.time);
                gam_o(k,:) = warp_f_gamma(gam_o(k,:),gamI_o,obj.time);
            end
            
            %% Aligned data & stats
            obj.fn = f1(:,:,r+1);
            obj.qn = q(:,:,r+1);
            obj.q0 = q(:,:,1);
            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            obj.mqn = mq(:,r+1);
            obj.fmean = mean(obj.f(1,:))+cumtrapz(t,obj.mqn.*abs(obj.mqn));
            
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean,gam_o(ii,:),t);
            end
            var_fgam = var(fgam,[],2);
            
            obj.stats.orig_var = trapz(t,std_f0.^2);
            obj.stats.amp_var = trapz(t,std_fn.^2);
            obj.stats.phase_var = trapz(t,var_fgam);
            
            obj.gam = gam_o.';
            [~,fy] = gradient(obj.gam,binsize,binsize);
            obj.psi = sqrt(fy+eps);
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            
            obj.qun = qun_o(1:r-1);
            obj.lambda = lambda;
            obj.method = option.method;
            obj.gamI = gamI_o;
            obj.rsamps = false;
            obj.type = 'median';
        end
        
        function obj = multiple_align_functions(obj, mu, lambda, option)
            % MULTIPLE_ALIGN_FUNCTIONS Group-wise function alignment to specified mean
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the elastic square-root
            % slope (srsf) framework.
            %
            % Usage:  obj.multiple_align_functions(mu)
            %         obj.multiple_align_functions(lambda)
            %         obj.multiple_align_functions(lambda, option)
            %
            % Input:
            % mu: vector of function to align to
            % lambda: regularization parameter
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.showplot = 1; % turns on and off plotting
            % option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS,expBayes)
            % option.w = 0.0; % BFGS weight
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output: structure containing
            % fdawarp object
            if nargin < 3
                lambda = 0;
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 1;
                option.method = 'DP1';
                option.w = 0.0;
                option.MaxItr = 20;
            elseif nargin < 4
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 1;
                option.method = 'DP1';
                option.w = 0.0;
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
            
            fprintf('\n lambda = %5.1f \n', lambda);
            
            [M, N] = size(obj.f);
            
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            
            
            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time);
            
            %% Compute the q-function of the plot
            mq = f_to_srvf(mu,obj.time);
            
            fn1 = zeros(M,N);
            qn1 = zeros(M,N);
            gam1 = zeros(N,size(q,1));
            mcmcopts.iter = option.MaxItr;
            mcmcopts.burnin = min(5e3,mcmcopts.iter/2);
            mcmcopts.alpha0 = 0.1;
            mcmcopts.beta0 = 0.1;
            tmp.betas = [0.5,0.5,0.005,0.0001];
            tmp.probs = [0.1,0.1,0.7,0.1];
            mcmcopts.zpcn = tmp;
            mcmcopts.propvar = 1;
            mcmcopts.initcoef = repelem(0, 20).';
            mcmcopts.npoints = 200;
            mcmcopts.extrainfo = true;
            if option.parallel == 1
                parfor k = 1:N
                    if (strcmpi(option.method,'expBayes'))
                        out_e = pairwise_align_bayes(mu, obj.f(:,k), obj.time, mcmcopts);
                        gam1(k,:) = out_e.gamma;
                    else
                        gam1(k,:) = optimum_reparam(mq,q(:,k),obj.time,lambda,option.method,option.w, ...
                            mu(1), obj.f(1,k));
                    end
                    fn1(:,k) = warp_f_gamma(obj.f(:,k,1),gam1(k,:),obj.time);
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time);
                end
            else
                for k = 1:N
                    if (strcmpi(option.method,'expBayes'))
                        out_e = pairwise_align_bayes(mu, obj.f(:,k), obj.time, mcmcopts);
                        gam1(k,:) = out_e.gamma;
                    else
                        gam1(k,:) = optimum_reparam(mq,q(:,k),obj.time,lambda,option.method,option.w, ...
                            mu(1), obj.f(1,k));
                    end
                    fn1(:,k) = warp_f_gamma(obj.f(:,k,1),gam1(k,:),obj.time);
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time);
                end
            end
            
            obj.gamI = SqrtMeanInverse(gam1);
            
            %% Aligned data & stats
            obj.q0 = q;
            obj.fn = fn1;
            obj.qn = qn1;
            obj.gam = gam1;
            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            obj.mqn = mq;
            obj.fmean = mean(obj.f(1,:))+cumtrapz(obj.time,obj.mqn.*abs(obj.mqn));
            
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean,obj.gam(ii,:),obj.time);
            end
            var_fgam = var(fgam,[],2);
            
            obj.stats.orig_var = trapz(obj.time,std_f0.^2);
            obj.stats.amp_var = trapz(obj.time,std_fn.^2);
            obj.stats.phase_var = trapz(obj.time,var_fgam);
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            
            obj.psi = [];
            obj.lambda = lambda;
            obj.method = option.method;
            obj.rsamps = false;
        end
        
        function obj = gauss_model(obj, n, sort_samples)
            % GAUSS_MODEL Gaussian gnerative model
            % -------------------------------------------------------------------------
            % This function models the functional data using a Gaussian model extracted
            % from the principal components of the srvfs
            %
            % Usage: obj.gauss_model(n, sort_samples)
            %
            % Inputs:
            % n: number of random samples (n = 1)
            % sort_samples: sort samples (e.g., = false)
            %
            % Output:
            % fdawarp object
            
            if (isempty(obj.type))
                error('Please align first');
            end
            if nargin < 2
                n = 1;
                sort_samples = false;
            elseif nargin < 3
                sort_samples = false;
            end
            %% Separated and Warped Data
            % sampling from the estimated model
            [M, ~] = size(obj.fn);
            binsize = mean(diff(obj.time));
            
            % compute mean and covariance in q-domain
            id = round(length(obj.time)/2);
            q_new = obj.qn;
            mq_new = mean(obj.qn,2);
            m_new = sign(obj.fn(id,:)).*sqrt(abs(obj.fn(id,:)));
            mqn2 = [mq_new; mean(m_new)];
            C = cov([q_new;m_new]');
            
            q_s = mvnrnd(mqn2', C, n)';
            
            % compute the correspondence to the original function domain
            f_s = zeros(M,n);
            for k = 1:n
                f_s(:,k) = cumtrapzmid(obj.time,q_s(1:(end-1),k).*abs(q_s(1:(end-1),k)),sign(q_s(end,k))*(q_s(end,k)^2),id);
            end
            fbar = mean(obj.fn,2);
            fsbar = mean(f_s,2);
            err = repmat(fbar-fsbar,1,n);
            f_s = f_s + err;
            
            % random warping generation
            rgam = randomGamma(obj.gam.',n);
            
            % combine samples
            if sort_samples
                %%%% sort functions and warpings
                mx = max(f_s);
                [~, seq1] = sort(mx);
                
                % compute the psi-function
                psi1 = zeros(n,size(rgam,2));
                len = zeros(1,n);
                ip = zeors(1,n);
                for i = 1:n
                    psi1(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
                    ip(i) = ones(1,M)*psi1(i,:)'/M;
                    len(i) = acos(ones(1,M)*psi1(i,:)'/M);
                end
                [~, seq2] = sort(len);
                
                % combine x-variability and y-variability
                f_c = zeros(size(obj.fn,1),n);
                for k = 1:n
                    f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,seq1(k)), invertGamma(rgam(seq2(k),:)')');
                    while sum(isnan(f_c(:,k))) >= 1
                        rgam2 = randomGamma(obj.gam.',1);
                        f_c(:,k) = warp_f_gamma(f_s(:,k), invertGamma(rgam2.'), (0:M-1)/(M-1));
                    end
                end
            else
                % compute the psi-function
                psi1 = zeros(n,size(rgam,2));
                len = zeros(1,n);
                ip = zeros(1,n);
                for i = 1:n
                    psi1(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
                    ip(i) = ones(1,M)*psi1(i,:)'/M;
                    len(i) = acos(ones(1,M)*psi1(i,:)'/M);
                end
                
                % combine x-variability and y-variability
                f_c = zeros(size(obj.fn,1),n);
                for k = 1:n
                    f_c(:,k) = interp1((0:M-1)/(M-1), f_s(:,k), invertGamma(rgam(k,:)')');
                    while sum(isnan(f_c(:,k))) >= 1
                        rgam2 = randomGamma(obj.gam.',1);
                        f_c(:,k) = warp_f_gamma(f_s(:,k), invertGamma(rgam2.'), (0:M-1)/(M-1));
                    end
                end
            end
            
            obj.fs = f_s;
            obj.gams = rgam;
            obj.ft = f_c;
            obj.qs = q_s(1:M,:);
            obj.rsamps = true;
        end
        
        function plot(obj)
            % plot plot functional alignment results
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            
            if (obj.rsamps)
                
                figure(1); clf
                M = length(obj.time);
                plot(obj.time, obj.ft, 'linewidth', 1);
                title('Random Functions', 'fontsize', 16);
                
                figure(2); clf;
                plot((0:M-1)/(M-1), obj.gams, 'linewidth', 1);
                axis square;
                title('Random Warping functions', 'fontsize', 16);
                
                figure(3); clf;
                plot(obj.time, obj.fs, 'LineWidth',1);
                title('Random Aligned Functions', 'fontsize', 16);
                
            else
                
                figure(1); clf;
                plot(obj.time, obj.f, 'linewidth', 1);
                title('Original data', 'fontsize', 16);
                
                if (~isempty(obj.gam))
                    mean_f0 = mean(obj.f, 2);
                    std_f0 = std(obj.f, 0, 2);
                    mean_fn = mean(obj.fn, 2);
                    std_fn = std(obj.fn, 0, 2);
                    figure(2); clf;
                    M = length(obj.time);
                    plot((0:M-1)/(M-1), obj.gam, 'linewidth', 1);
                    axis square;
                    title('Warping functions', 'fontsize', 16);
                    
                    figure(3); clf;
                    plot(obj.time, obj.fn, 'LineWidth',1);
                    title(['Warped data, \lambda = ' num2str(obj.lambda)], 'fontsize', 16);
                    
                    figure(4); clf;
                    plot(obj.time, mean_f0, 'b-', 'linewidth', 1); hold on;
                    plot(obj.time, mean_f0+std_f0, 'r-', 'linewidth', 1);
                    plot(obj.time, mean_f0-std_f0, 'g-', 'linewidth', 1);
                    title('Original data: Mean \pm STD', 'fontsize', 16);
                    
                    figure(5); clf;
                    plot(obj.time, mean_fn, 'b-', 'linewidth', 1); hold on;
                    plot(obj.time, mean_fn+std_fn, 'r-', 'linewidth', 1);
                    plot(obj.time, mean_fn-std_fn, 'g-', 'linewidth', 1);
                    title(['Warped data, \lambda = ' num2str(obj.lambda) ': Mean \pm STD'], 'fontsize', 16);
                    
                    figure(6); clf;
                    plot(obj.time, obj.fmean, 'g','LineWidth',1);
                    title(['f_{mean}, \lambda = ' num2str(obj.lambda)], 'fontsize', 16);
                end
            end
        end
        
    end
    
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
