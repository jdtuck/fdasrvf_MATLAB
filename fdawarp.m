classdef fdawarp
    % fdawarp elastic fda functional class
    %   fdawarp object contains the ability to align and plot functional
    %   data and is required for follow on analysis
    
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
        type   % alignment type
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
            % Usage:  out = time_warping(f,t)
            %         out = time_warping(f,t,lambda)
            %         out = time_warping(f,t,lambda,option)
            %
            % Input:
            % f (M,N): matrix defining N functions of M samples
            % t : time vector of length M
            % lambda: regularization parameter
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.showplot = 1; % turns on and off plotting
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
            % Input:
            % f (M,N): matrix defining N functions of M samples
            % t : time vector of length M
            % lambda: regularization parameter
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.showplot = 1; % turns on and off plotting
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
                option.showplot = 1;
                option.method = 'DP1';
                option.w = 0.0;
                option.MaxItr = 20;
            elseif nargin < 3
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
            %% Parameters
            
            fprintf('\n lambda = %5.1f \n', lambda);
            
            t = obj.time;
            
            binsize = mean(diff(t));
            [M, N] = size(obj.f);
            f1 = obj.f;
            
            if option.smooth == 1
                f1 = smooth_data(f1, option.sparam);
            end
            
            if option.showplot == 1
                figure(1); clf;
                plot(t, obj.f, 'linewidth', 1);
                title('Original data', 'fontsize', 16);
                pause(0.1);
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
            mq(:,r+1) = warp_q_gamma(mq(:,r),gamI_o,t);
            for k = 1:N
                q(:,k,r+1) = warp_q_gamma(q(:,k,r),gamI_o,t);
                f1(:,k,r+1) = warp_q_gamma(f1(:,k,r),gamI_o,t);
                gamI_o(k,:) = interp1(t, gamI_o(k,:), (t(end)-t(1)).*gamI_o + t(1));
            end
            
            %% Aligned data & stats
            obj.fn = f1(:,:,r+1);
            obj.qn = q(:,:,r+1);
            obj.q0 = q(:,:,1);
            std_f0 = std(f0, 0, 2);
            std_fn = std(fn, 0, 2);
            obj.mqn = mq(:,r+1);
            fmedian = mean(obj.f(1,:))+cumtrapz(t,obj.mqn.*abs(obj.mqn));
            
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(fmedian,gam_o(ii,:),t);
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
        
        function out = Amplitude_Boxplot(out_warp, alpha, k_a, figs)
            % AMPLITUDE_BOXPLOT Functional Amplitude Boxplot
            % -------------------------------------------------------------------------
            %
            % This function constructs the amplitude boxplot
            %
            % Usage:  out = Amplitude_Boxplot(out_warp, alpha, k_a, figs)
            %
            % Input:
            % out_warp: struct from time_warping_median of aligned data using the median
            % alpha: quantile value (e.g.,=.05, i.e., 95\%)
            % ka: scalar for outlier cutoff (e.g.,=1)
            % figs: shows plots of functions (e.g., = true)
            %
            % Output: structure containing
            % median_y: median function
            % Q1: First quartile
            % Q3: Second quartile
            % Q1a: First quantile based on alpha
            % Q3a: Second quantile based on alpha
            % minn: minimum extreme function
            % maxx: maximum extreme function
            % outlier_index: indexes of outlier functions
            % fmedian: median function
            
            if (~strcmpi(obj.type,'median'))
                error('rerun alignment with time_warping_median');
            end
            f_tilde = obj.fn;
            f_median = obj.fmedian;
            q_tilde = obj.qn;
            q_median = obj.mqn;
            t = obj.time;
            
            [M, N] = size(f_tilde);
            lambda1 = 0.5;
            
            % amplitude median
            median_y = f_median;
            
            % compute amplitude distances
            dy = zeros(1,N);
            for i = 1:N
                dy(i) = sqrt(trapz(t,(q_median-q_tilde(:,i)).^2));
            end
            [~, dy_ordering] = sort(dy);
            CR_50 = dy_ordering(1:ceil(N/2));       % 50% Central Region
            m = max(dy(CR_50));                     % Maximal amplitude distance within 50% Central Region
            
            % identify amplitude quartiles
            angle = zeros(length(CR_50), length(CR_50));
            energy = zeros(length(CR_50), length(CR_50));
            for i = 1:(length(CR_50)-1)
                for j = (i+1):length(CR_50)
                    q1 = q_tilde(:,CR_50(i)) - q_median;
                    q3 = q_tilde(:,CR_50(j)) - q_median;
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda1) * (dy(CR_50(i))/m + dy(CR_50(j))/m) - lambda1 * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            Q1_index = CR_50(maxloc_row);
            Q3_index = CR_50(maxloc_col);
            Q1_q = q_tilde(:,Q1_index);
            Q3_q = q_tilde(:,Q3_index);
            Q1 = f_tilde(:,Q1_index);
            Q3 = f_tilde(:,Q3_index);
            
            % identify amplitude quantiles
            [~, dy_ordering] = sort(dy);
            CR_alpha = dy_ordering(1:round(N*(1-alpha)));       % (1-alpha)% Central Region
            m = max(dy(CR_alpha));                     % Maximal amplitude distance within (1-alpha)% Central Region
            angle = zeros(length(CR_alpha), length(CR_alpha));
            energy = zeros(length(CR_alpha), length(CR_alpha));
            for i = 1:(length(CR_alpha)-1)
                for j = (i+1):length(CR_alpha)
                    q1 = q_tilde(:,CR_alpha(i)) - q_median;
                    q3 = q_tilde(:,CR_alpha(j)) - q_median;
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda1) * (dy(CR_alpha(i))/m + dy(CR_alpha(j))/m) - lambda1 * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            Q1a_index = CR_alpha(maxloc_row);
            Q3a_index = CR_alpha(maxloc_col);
            Q1a_q = q_tilde(:,Q1a_index);
            Q3a_q = q_tilde(:,Q3a_index);
            Q1a = f_tilde(:,Q1a_index);
            Q3a = f_tilde(:,Q3a_index);
            
            % compute amplitude whiskers
            IQR = dy(Q1_index)+dy(Q3_index);
            v1 = Q1_q - q_median;
            v3 = Q3_q - q_median;
            upper_q = Q3_q + k_a * IQR * v3 / sqrt(trapz(t,v3.^2));
            lower_q = Q1_q + k_a * IQR * v1 / sqrt(trapz(t,v1.^2));
            
            upper_dis = sqrt(trapz(t,(upper_q - q_median).^2));
            lower_dis = sqrt(trapz(t,(lower_q - q_median).^2));
            whisker_dis = max([lower_dis upper_dis]);
            
            % identify amplitude outliers
            outlier_index = [];
            for i = 1:N
                if dy(dy_ordering(N+1-i)) > whisker_dis
                    outlier_index = [outlier_index; dy_ordering(N+1-i)];
                else
                    break
                end
            end
            
            % identify amplitude extremes
            distance_to_upper=inf(1,N);
            distance_to_lower=inf(1,N);
            out_50_CR = setdiff(setdiff((1:N), CR_50), outlier_index);
            for i = 1:length(out_50_CR)
                j = out_50_CR(i);
                distance_to_upper(j) = sqrt(trapz(t,(upper_q - q_tilde(:,j)).^2));
                distance_to_lower(j) = sqrt(trapz(t,(lower_q - q_tilde(:,j)).^2));
            end
            [~, max_index] = min(distance_to_upper);
            [~, min_index] = min(distance_to_lower);
            min_q = q_tilde(:,min_index);
            max_q = q_tilde(:,max_index);
            minn = f_tilde(:,min_index);
            maxx = f_tilde(:,max_index);
            
            s = linspace(0,1,100);
            Fs2 = zeros(length(t), 595);
            Fs2(:,1) = (1-s(1)) * minn + s(1) * Q1;    % Final surface plot
            for j=2:100
                Fs2(:,j) = (1-s(j)) * minn + s(j) * Q1a;
                Fs2(:,99+j) = (1-s(j)) * Q1a + s(j) * Q1;
                Fs2(:,198+j) = (1-s(j)) * Q1 + s(j) * f_median;
                Fs2(:,297+j) = (1-s(j)) * f_median + s(j) * Q3;
                Fs2(:,396+j) = (1-s(j)) * Q3 + s(j) * Q3a;
                Fs2(:,495+j) = (1-s(j)) * Q3a + s(j) * maxx;
            end
            d1=sqrt(trapz(t,(q_median-Q1_q).^2));
            d1a=sqrt(trapz(t,(Q1_q-Q1a_q).^2));
            dl=sqrt(trapz(t,(Q1a_q-min_q).^2));
            d3=sqrt(trapz(t,(q_median-Q3_q).^2));
            d3a=sqrt(trapz(t,(Q3_q-Q3a_q).^2));
            du=sqrt(trapz(t,(Q3a_q-max_q).^2));
            part1=linspace(-d1-d1a-dl,-d1-d1a,100);
            part2=linspace(-d1-d1a,-d1,100);
            part3=linspace(-d1,0,100);
            part4=linspace(0,d3,100);
            part5=linspace(d3,d3+d3a,100);
            part6=linspace(d3+d3a,d3+d3a+du,100);
            allparts=[part1,part2(2:100),part3(2:100),part4(2:100),part5(2:100),part6(2:100)];
            [U,V]=meshgrid(t,allparts);
            U=U';
            V=V';
            
            plt.U=U;
            plt.V=V;
            plt.Fs2 = Fs2;
            plt.allparts = allparts;
            plt.d1 = d1;
            plt.d1a = d1a;
            plt.dl = dl;
            plt.d3 = d3;
            plt.d3a = d3a;
            plt.du = du;
            plt.Q1q = Q1a_q;
            plt.Q3q = Q3a_q;
            
            if (figs)
                figure(310); clf;
                plot(t, f_median, 'black','linewidth', 2);
                hold on;
                plot(t, Q1, 'blue','linewidth', 2);
                plot(t, Q3, 'blue', 'linewidth', 2);
                plot(t, Q1a, 'green', 'linewidth', 2);
                plot(t, Q3a, 'green', 'linewidth', 2);
                plot(t,minn,'red', 'linewidth',2);
                plot(t,maxx,'red', 'linewidth',2);
                xlim([t(1) t(end)]);
                ylim auto;
                
                
                figure(311); clf;
                surf(U,V,Fs2);
                hold on;
                shading flat;
                plot3(t,zeros(1,M),f_median,'k','LineWidth',3)
                plot3(t,repmat(-d1,M,1),Q1,'b','LineWidth',3)
                plot3(t,repmat(-d1-d1a,M,1),Q1a,'g','LineWidth',3)
                plot3(t,repmat(-d1-d1a-dl,M,1),minn,'r','LineWidth',3)
                plot3(t,repmat(d3,M,1),Q3,'b','LineWidth',3)
                plot3(t,repmat(d3+d3a,M,1),Q3a,'g','LineWidth',3)
                plot3(t,repmat(d3+d3a+du,M,1),maxx,'r','LineWidth',3)
            end
            
            out.f_median = median_y;
            out.q_median = q_median;
            out.Q1 = Q1;
            out.Q3 = Q3;
            out.Q1a = Q1a;
            out.Q3a = Q3a;
            out.minn = minn;
            out.maxx = maxx;
            out.outlier_index = outlier_index;
            out.plt = plt;
        end
        
        function plot(obj)
            % plot plot functional alignment results
            % -------------------------------------------------------------------------
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
