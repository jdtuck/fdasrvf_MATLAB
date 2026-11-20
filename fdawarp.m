classdef fdawarp
    %fdawarp A class to provide alignment of functional data using SRVF
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
    %    lambda - warping penalty
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
    %   time_warping_bayes - align functions and find karcher mean using Bayesian
    %   approach
    %   time_warping_bayes_infHMC - Bayesian group alignment using infinite-dimensional HMC
    %   time_warping_median - align functions and find karcher median
    %   multiple_align_functions - align functions to a specified mean
    %   gauss_model - model functions using a generative gaussian model on
    %       fPCA space
    %   align_new - align new functions to the fitted template
    %   plot - plot results and functions in object
    %   plot_mcmc - plot MCMC diagnostics after time_warping_bayes
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
        penalty % warping penalty
        gamI   % invserse warping function
        rsamps % random samples
        fs     % random aligned functions
        gams   % random warping functions
        ft     % random warped functions
        qs     % random aligned srvfs
        type   % alignment type
        mcmc   % mcmc output if bayesian
        lam_opt % optimum lambda from ppd

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

            if (size(f,1) ~= length(time))
                error('Columns of f and time must be equal')
            end

            obj.f = f;
            obj.time = time;
        end

        function obj = ppd(obj, max_lam, num_lam, pt, srvf, mu, option)
            % PPD Compute Peak Persistance Diagram
            % -------------------------------------------------------------------------
            % This computes the peak persistance diagram over a range of
            % lambda. This can help determine the proper elasticity (penalty).
            % This can be slow and recommended to run in parallel
            %
            % Usage:  obj.ppd()
            %         obj.ppd(max_lam, num_lam)
            %         obj.ppd(max_lam, num_lam, pt)
            %         obj.ppd(max_lam, num_lam, pt, option)
            %
            % Input:
            % max_lam: max regularization parameter (default 2)
            % num_lam: number of steps (default 10)
            % pt: the percentile of negative curvature of raw data (default
            % .15)
            % srvf: perfrom in srvf space (default F)
            % mu: provided mean to align (default NaN)
            %
            % default options
            % option.parallel = 0; % turns on MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
            % option.penalty = 'roughness'; % warping penalty weighted by lambda (none, roughness, l2gam, l2psi, geodesic)
            % option.spl = true; % use interpolating spline derivative (false: smoothing spline, lossy)
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdawarp object
            arguments
                obj
                max_lam = 2
                num_lam = 10
                pt = 0.15
                srvf = false
                mu = NaN
                option.parallel = 1;
                option.nworkers = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP';
                option.penalty = 'roughness';
                option.spl = true;
                option.MaxItr = 20;
            end

            % Obtain a lambda candidate set
            lam_vec = linspace(0, max_lam, num_lam);
            fns = cell(num_lam, 1);
            for i = 1:num_lam
                if isnan(mu)
                    obj = obj.time_warping(lam_vec(i), option);
                else
                    obj = obj.multiple_align_functions(mu, lam_vec(i), ...
                        parallel=option.parallel,nworkers=option.nworkers,closepool=option.closepool,...
                        smooth=option.smooth,sparam=option.sparam,...
                        method=option.method,penalty=option.penalty,spl=option.spl,MaxItr=option.MaxItr);
                end
                if srvf
                    fns{i} = obj.qn;
                else
                    fns{i} = obj.fn;
                end
            end
            if srvf
                f0 = obj.q0;
            else
                f0 = obj.f;
            end

            % Peak Persistent Diagrams
            % Get the threshold for significant peak
            diff_t = mean(diff(obj.time));
            taus = [];
            
            % Compute tau values
            for i = 1:size(obj.f, 2)
                idx = islocalmax(f0(:,i));
                df2 = gradient(gradient(f0(:,i), diff_t), diff_t);
                tau = -df2 / max(-df2);
                tau = max(tau,0);
                taus = [taus, tau(idx)'];
            end

            th = prctile(taus, pt * 100);

            [IndicatorMatrix, ~, Heights, Locs, Labels, FNm] = getPPDinfo(obj.time, fns, lam_vec, th);
            persistent_peak_labels = getPersistentPeaks(IndicatorMatrix');
                
            % Choose optimal lambda
            ref_row = zeros(1,size(IndicatorMatrix,2));
            ref_row(persistent_peak_labels) = 1;
            
            n_lams = size(IndicatorMatrix, 1);
            hamming_distances = zeros(1, n_lams);
            exact_match_indices = [];

            for i = 1:n_lams
                comp = double(~isnan(IndicatorMatrix(i, :)));
                
                % Check for exact match
                if isequal(comp, ref_row)
                    exact_match_indices = [exact_match_indices, i];
                end    
                % Calculate Hamming distance
                hamming_distances(i) = sum(comp ~= ref_row);
            end
            
            if ~isempty(exact_match_indices)
                idx_opt = min(exact_match_indices);
            else
                min_distance = min(hamming_distances);
                
                min_distance_indices = find(hamming_distances == min_distance);
                idx_opt = min(min_distance_indices);
            end

            % Draw PPD BarChart
            drawPPDBarChart(IndicatorMatrix,Heights,lam_vec,idx_opt)
            % Draw PPD Surface
            drawPPDSurface(obj.time,lam_vec,FNm,Heights,Locs,IndicatorMatrix,Labels,idx_opt)
            obj.lam_opt = lam_vec(idx_opt);
            
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
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
            % option.penalty = 'roughness'; % warping penalty weighted by lambda (none, roughness, l2gam, l2psi, geodesic)
            % option.spl = true; % use interpolating spline derivative (false: smoothing spline, lossy)
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdawarp object
            arguments
                obj
                lambda = 0;
                option.parallel = 0;
                option.nworkers = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP';
                option.penalty = 'roughness';
                option.spl = true;
                option.MaxItr = 20;
            end

            % time warping on a set of functions
            if option.parallel == 1
                start_parpool(option.nworkers);
            end
            %% Parameters

            fprintf('\n lambda = %5.2f \n', lambda);

            binsize = mean(diff(obj.time));
            [M, N] = size(obj.f);

            f1 = obj.f;
            if option.smooth == 1
                f1 = smooth_data(f1, option.sparam);
            end

            %% Compute the q-function of the plot
            q = f_to_srvf(f1,obj.time,~option.spl,option.spl,option.parallel,option.nworkers);

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
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                        mf(1), f1(1,k,1), 10, option.penalty);
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                        mf(1), f1(1,k,1), 10, option.penalty);
                end
            end
            gamI_o = SqrtMeanInverse(gam_o');
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
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                            mf(1,r), f1(1,k,1), 10, option.penalty);
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time,~option.spl,option.spl);
                    end
                else
                    for k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                            mf(1,r), f1(1,k,1), 10, option.penalty);
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time,~option.spl,option.spl);
                    end
                end
                q(:,:,r+1) = q_temp;
                f1(:,:,r+1) = f_temp;

                ds(r+1) = sum(simps(obj.time, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
                    lambda*sum(simps(obj.time, (1-sqrt(max(gam_dev',0))).^2));

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
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                        mf(1,r), f1(1,k,1), 10, option.penalty);
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method, ...
                        mf(1,r), f1(1,k,1), 10, option.penalty);
                end
            end
            gamI_o = SqrtMeanInverse(gam_o');
            % Center gam, then warp the original data once by the final
            % gam, so fn = f o gam and mqn = mean(qn) hold exactly.
            % (Previously fn/qn kept the previous iteration's alignment
            % while gam came from this final matching step.)
            for k = 1:N
                gam_o(k,:) = interp1(obj.time, gam_o(k,:), (obj.time(end)-obj.time(1)).*gamI_o + obj.time(1));
                f1(:,k,r+1) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                q(:,k,r+1) = f_to_srvf(f1(:,k,r+1),obj.time,~option.spl,option.spl);
            end
            mq(:,r+1) = mean(q(:,:,r+1),2);

            %% Aligned data & stats
            obj.fn = f1(:,:,r+1);
            obj.qn = q(:,:,r+1);
            obj.q0 = q(:,:,1);
            std_f0 = std(f1(:,:,1), 0, 2);
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
            obj.penalty = option.penalty;
            obj.gamI = gamI_o;
            obj.rsamps = false;
            obj.mcmc = [];
            obj.type = 'mean';
        end

        function obj = time_warping_bayes(obj, mcmcopts)
            % TIME_WARPING_BAYES Align multiple functions using Bayesian method
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the elastic square-root
            % slope (srsf) framework. This method uses a Bayesian approach.
            % It is based on mapping warping functions to a hypersphere, and a
            % subsequent exponential mapping to a tangent space. In the tangent space,
            % the Z-mixture pCN algorithm is used to explore both local and global
            % structure in the posterior distribution.
            %
            % The Z-mixture pCN algorithm uses a mixture distribution for the proposal
            % distribution, controlled by input parameter zpcn. The zpcn$betas must be
            % between 0 and 1, and are the coefficients of the mixture components, with
            % larger coefficients corresponding to larger shifts in parameter space. The
            % zpcn$probs give the probability of each shift size.
            %
            % Usage:  obj.time_warping_bayes()
            %         obj.time_warping_bayes(mcmcopts)
            %
            % Input:
            %
            % default mcmc options
            % mcmcopts.iter = 2000; % number of iterations
            % mcmcopts.burnin = min(5e3,mcmcopts.iter/2); % number of burnin
            % mcmcopts.alpha0 = 0.1; # inverse gamma prior parameters
            % mcmcopts.beta0 = 0.1; # inverse gamma prior parameters
            % tmp.betas = [0.5,0.5,0.005,0.0001]; %pczn paramaters
            % tmp.probs = [0.1,0.1,0.7,0.1];
            % mcmcopts.zpcn = tmp;
            % mcmcopts.propvar = 1; % proposal variance
            % mcmcopts.initcoef = zeros(20, size(f,2)) % init warping function coef
            % mcmcopts.nbasis = 40; % number of basis elements for q
            % mcmcopts.npoints = 200; % number of sample interpolation points
            % mcmcopts.extrainfo = true; % return extra info about mcmc
            % mcmcopts.ncenter = 100; % number of iterations inbetween centering gam
            %
            % Output:
            % fdawarp object
            %
            % if extrainfo object contains
            % mcmc.accept: accept of q samples
            % mcmc.betas_ind
            % mcmc.gamma_mat: posterior gammas
            % mcmc.qstar_mat: posterior q_star
            % mcmc.fmean_mat: posterior fmean
            % mcmc.gamma_stats: posterior gamma stats
            % mcmc.qstar_stats: posterior q_star stats
            % mcmc.fmean_stats: posterior fmean stats
            arguments
                obj
                mcmcopts.iter = 2000;
                mcmcopts.burnin = min(5e3,2000/2);
                mcmcopts.alpha0 = 0.1;
                mcmcopts.beta0 = 0.1;
                mcmcopts.zpcn = struct('betas',[0.5,0.5,0.005,0.0001],'probs',[0.1,0.1,0.7,0.1]);
                mcmcopts.propvar = 4;
                mcmcopts.initcoef = zeros(20, size(obj.f,2));
                mcmcopts.nbasis = 30;
                mcmcopts.npoints = 200;
                mcmcopts.extrainfo = true;
                mcmcopts.ncenter = 100;
            end
            [M, N] = size(obj.f);

            if (~iscolumn(obj.time))
                obj.time = obj.time.';
            end

            if (length(mcmcopts.zpcn.betas) ~= length(mcmcopts.zpcn.probs))
                error('In zpcn, betas must equal length of probs')
            end

            if (mod(size(mcmcopts.initcoef), 1) ~= 0)
                error('Length of mcmcopts.initcoef must be even')
            end

            if (mod(mcmcopts.nbasis, 1) ~= 0)
                error('mcmcopts.nbasis must be even')
            end

            % Number of sig figs to report in gamma_mat
            SIG_GAM = 13;
            iter = mcmcopts.iter;

            % normalize timet to [0,1]
            % ([a,b] - a) / (b-a) = [0,1]
            obj.time = (obj.time-min(obj.time))./(max(obj.time)-min(obj.time));

            % parameter settings
            pw_sim_global_burnin = mcmcopts.burnin;
            valid_index = pw_sim_global_burnin:iter;
            pw_sim_global_Mg = size(mcmcopts.initcoef,1)/2;
            g_coef_ini = mcmcopts.initcoef;
            numSimPoints = mcmcopts.npoints;
            pw_sim_global_domain_par = linspace(0,1,numSimPoints).';
            g_basis = basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mg, 1);
            pw_sim_global_Mq = mcmcopts.nbasis/2;
            q_basis = basis_bspline(pw_sim_global_domain_par, 2*pw_sim_global_Mq, 3); %basis_fourier(pw_sim_global_domain_par, pw_sim_global_Mq, 2);
            sigma1_ini = 1;
            zpcn = mcmcopts.zpcn;

            % srsf transformation
            f0 = zeros(length(g_basis.x),N);
            for ii = 1:N
                f0(:,ii) = interp1(obj.time,obj.f(:,ii),g_basis.x,'linear');
            end
            qo = f_to_srvf(f0, g_basis.x);

            pw_sim_global_sigma_g = mcmcopts.propvar;
            pw_sim_global_sigma_q = (quantile(qo(:),.95)-quantile(qo(:),.05))*5;
            pw_sim_global_sigma_q_const = (quantile(qo(:),.95)-quantile(qo(:),.05))*.5;

            function result = propose_g_coef(g_coef_curr)
                pCN_beta = zpcn.betas;
                pCN_prob = zpcn.probs;
                probm = [0, cumsum(pCN_prob)];
                z = rand;
                for i = 1:length(pCN_beta)
                    if (z <= probm(i+1) && z > probm(i))
                        g_coef_new = normrnd(0, pw_sim_global_sigma_g ./ repelem(1:pw_sim_global_Mg,2), 1, pw_sim_global_Mg * 2);
                        result.prop = sqrt(1-pCN_beta(i)^2) * g_coef_curr + pCN_beta(i) * g_coef_new.';
                        result.ind = i;
                    end
                end
            end

            function result = propose_q_star(q_star_curr)
                pCN_beta = zpcn.betas;
                pCN_prob = zpcn.probs;
                probm = [0, cumsum(pCN_prob)];
                z = rand;
                sdvec = pw_sim_global_sigma_q./repelem(1:pw_sim_global_Mq,2);
                for i = 1:length(pCN_beta)
                    if (z <= probm(i+1) && z > probm(i))
                        q_star_new = normrnd(0, sdvec, 1, pw_sim_global_Mq * 2);
                        result.prop = sqrt(1-pCN_beta(i)^2) * q_star_curr + pCN_beta(i) * q_star_new.';
                        result.ind = i;
                    end
                end
            end

            for ii = 1:N
                tmp = f_exp1(f_basistofunction(g_basis.x,0,g_coef_ini(:,ii),g_basis, false));
                if (min(tmp.y)<0)
                    error("Invalid initial value of g for sample %d", ii)
                end
            end

            % result objects
            result.g_coef = zeros(iter,size(g_coef_ini,1),size(g_coef_ini,2));
            result.q_star = zeros(iter,length(g_basis.x));
            result.sigma1 = zeros(1,iter);
            result.logl = zeros(1,iter);
            result.SSE = zeros(iter,N);
            result.accept = zeros(1,iter);
            result.accept_betas = zeros(1,iter);

            % init
            g_coef_curr = g_coef_ini;
            sigma1_curr = sigma1_ini;
            mnq = mean(qo,2);
            dqq = sqrt(sum((qo - mnq*ones(1,N)).^2,1));
            [~, min_ind] = min(dqq);
            q_star_curr = qo(:,min_ind);
            q.y = qo;
            q.x = g_basis.x;
            q_star_ini = q_star_curr;
            SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false), ...
                q,q_star_ini);
            logl_curr = f_logl(f_basistofunction(g_basis.x,0,g_coef_ini,g_basis,false), ...
                q,q_star_ini,sigma1_ini^2,SSE_curr);

            result.g_coef(1,:,:) = g_coef_ini;
            result.sigma1(1) = sigma1_ini;
            result.q_star(1,:) = q_star_curr;
            result.SSE(1,:) = SSE_curr;
            result.logl(1) = logl_curr;

            q_star_coef_curr = q_basis.matrix\q_star_curr;
            clear q_star_curr;
            q_star_curr.x = g_basis.x;
            q_star_curr.y = q_star_ini;

            % update the chain for iter-1 times
            barobj = ProgressBar(iter-1, ...
                'Title', 'Running MCMC' ...
                );
            for m = 2:iter
                % update g
                for ii = 1:N
                    q_tmp = q;
                    q_tmp.y = q_tmp.y(:,ii);
                    [g_coef_curr(:,ii), ~, ~] = f_updateg_pw(g_coef_curr(:,ii), g_basis, sigma1_curr^2, q_tmp, q_star_curr.y, @propose_g_coef);
                end
                SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
                    q,q_star_curr.y);
                % update q_star
                [q_star_coef_curr, accept, zpcnInd] = f_updateq(g_coef_curr, g_basis, sigma1_curr^2, q, q_star_coef_curr, q_basis, SSE_curr, @propose_q_star);
                q_star_curr = f_basistofunction(q_basis.x,0,q_star_coef_curr,q_basis,false);

                % update sigma1
                SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
                    q,q_star_curr.y);
                newshape = N*length(obj.time)/2 + mcmcopts.alpha0;
                newscale = 1/2 * sum(SSE_curr) + mcmcopts.beta0;
                sigma1_curr = sqrt(1/gamrnd(newshape,1/newscale));
                logl_curr = f_logl(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
                    q,q_star_curr.y,sigma1_curr^2,SSE_curr);
                
                % center
                if ((mcmcopts.ncenter ~=0) && (mod(m,mcmcopts.ncenter)==0))
                    result_posterior_gamma = zeros(length(g_basis.x),N);
                    for ii = 1:N
                        g_temp = f_basistofunction(pw_sim_global_domain_par, 0, g_coef_curr(:,ii), g_basis, false);
                        psi_temp = f_exp1(g_temp);
                        % transform posterior mean of psi to gamma
                        tmp_gamma = f_phiinv(psi_temp);
                        gam0 = tmp_gamma.y;
                        tmp_gamma.y = norm_gam(gam0);
                        result_posterior_gamma(:,ii) = norm_gam(gam0);
                    end
                    gamI_o = SqrtMeanInverse(result_posterior_gamma);
                    q_star_curr.y = warp_q_gamma(q_star_curr.y,gamI_o,q_star_curr.x);
                    for k = 1:N
                        result_posterior_gamma(:,k) = interp1(q_star_curr.x, result_posterior_gamma(:,k), (obj.time(end)-obj.time(1)).*gamI_o + obj.time(1));
                        gamma.x = g_basis.x;
                        gamma.y = result_posterior_gamma(:,k);
                        psi1 = f_phi(gamma);
                        [~,yy] = f_exp1inv(psi1);
                        g_coef_curr(:,k) = g_basis.matrix\yy;
                    end
                    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
                    q,q_star_curr.y);
                    logl_curr = f_logl(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), ...
                    q,q_star_curr.y,sigma1_curr^2,SSE_curr);
                end

                % save update to results
                result.g_coef(m,:,:) = g_coef_curr;
                result.q_star(m,:) = q_star_curr.y;
                result.sigma1(m) = sigma1_curr;
                result.SSE(m,:) = SSE_curr;
                if (mcmcopts.extrainfo)
                    result.logl(m) = logl_curr;
                    result.accept(m) = accept;
                    result.accept_betas(m) = zpcnInd;
                end
                barobj.step([], [], []);
            end
            barobj.release();

            % calculate posterior mean of psi
            pw_sim_est_psi_matrix = zeros(length(pw_sim_global_domain_par), length(valid_index),N);
            for k = 1:length(valid_index)
                for ii = 1:N
                    g_temp = f_basistofunction(pw_sim_global_domain_par, 0, result.g_coef(valid_index(k),:,ii).', g_basis, false);
                    psi_temp = f_exp1(g_temp);
                    pw_sim_est_psi_matrix(:,k,ii) = psi_temp.y;
                end
            end

            result_posterior_psi_simDomain = cell(N,1);
            result_posterior_psi = zeros(M,N);
            result_posterior_gamma = zeros(M,N);
            f_warped = zeros(M,N);
            for ii = 1:N
                result_posterior_psi_simDomain{ii} = f_psimean(pw_sim_global_domain_par, pw_sim_est_psi_matrix(:,:,ii));
                % resample to same number of points as the input f1 and f2
                result_i = interp1(result_posterior_psi_simDomain{ii}.x, result_posterior_psi_simDomain{ii}.y, obj.time, 'linear', 'extrap');
                tmp_psi.x=obj.time;
                tmp_psi.y=result_i;
                result_posterior_psi(:,ii) = tmp_psi.y;
                % transform posterior mean of psi to gamma
                tmp_gamma = f_phiinv(tmp_psi);
                gam0 = tmp_gamma.y;
                tmp_gamma.y = norm_gam(gam0);
                result_posterior_gamma(:,ii) = norm_gam(gam0);
                f_warped(:,ii) = warp_f_gamma(obj.f(:,ii), result_posterior_gamma(:,ii), tmp_gamma.x);
            end

            if (mcmcopts.extrainfo)
                % matrix of posterior draws from gamma
                gamma_mat = zeros(length(obj.time),size(pw_sim_est_psi_matrix,2),N);
                gamma_stats = zeros(2,size(pw_sim_est_psi_matrix,2),N);
                for jj = 1:N
                    for ii = 1:size(pw_sim_est_psi_matrix,2)
                        result_i = interp1(result_posterior_psi_simDomain{jj}.x, pw_sim_est_psi_matrix(:,ii,jj), obj.time, 'linear', 'extrap');
                        result_i2.y=result_i;
                        result_i2.x=obj.time;
                        tmp = f_phiinv(result_i2);
                        gamma_mat(:,ii,jj) = round(norm_gam(tmp.y),SIG_GAM);
                        gamma_stats(:,ii,jj) = statsFun(gamma_mat(:,ii,jj));
                    end
                end
            end

            % return object
            obj.fn = f_warped;
            obj.gam = result_posterior_gamma;
            obj.psi = result_posterior_psi;
            obj.qn = f_to_srvf(f_warped,obj.time);
            obj.q0 = f_to_srvf(obj.f,obj.time);
            mqn_temp = mean(result.q_star(valid_index,:),1)';
            obj.mqn = interp1(g_basis.x, mqn_temp, obj.time, 'linear', 'extrap'); %mean(obj.qn,2);
            obj.fmean = mean(obj.f(1,:))+cumtrapz(obj.time,obj.mqn.*abs(obj.mqn));
            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean,obj.gam(:,ii),obj.time);
            end
            var_fgam = var(fgam,[],2);

            obj.stats.orig_var = trapz(obj.time,std_f0.^2);
            obj.stats.amp_var = trapz(obj.time,std_fn.^2);
            obj.stats.phase_var = trapz(obj.time,var_fgam);

            obj.qun = result.logl;
            obj.method = "bayesian";
            obj.gamI = interp1(g_basis.x, gamI_o, obj.time, 'linear', 'extrap');
            obj.rsamps = false;
            obj.type = 'mean';

            obj.mcmc = struct();
            obj.mcmc.g_coef = result.g_coef;
            obj.mcmc.sigma1 = result.sigma1;
            obj.mcmc.logl = result.logl;
            obj.mcmc.burnin = pw_sim_global_burnin;
            obj.mcmc.iter = iter;
            obj.mcmc.extrainfo = mcmcopts.extrainfo;

            if (mcmcopts.extrainfo)
                obj.mcmc.accept = result.accept(2:end);
                obj.mcmc.betas_ind = result.accept_betas(2:end);
                obj.mcmc.gamma_mat = gamma_mat;
                obj.mcmc.gamma_stats = gamma_stats;
                obj.mcmc.q_star_mat = result.q_star(valid_index,:)';
                obj.mcmc.q_star_stats = quantile(obj.mcmc.q_star_mat,[.025,0.975],2);
                obj.mcmc.fmean_mat = srvf_to_f(obj.mcmc.q_star_mat,q_basis.x,repelem(mean(obj.f(1,:)),length(valid_index)));
                obj.mcmc.fmean_stats = quantile(obj.mcmc.fmean_mat,[.025,0.975],2);
                obj.mcmc.zpcn = zpcn;
            end

        end

        function obj = time_warping_bayes_infHMC(obj, mcmcopts)
            % TIME_WARPING_BAYES_INFHMC Bayesian group-wise alignment using
            % infinite-dimensional Hamiltonian Monte Carlo
            % -------------------------------------------------------------------------
            % Hierarchical Bayesian model in SRVF space. For i = 1, ..., N
            %
            %   q*(t_k) = (q_i o gam_i)(t_k) sqrt(gam_i'(t_k)) + e_ik,  e_ik ~ N(0, sigma^2)
            %
            % with template q* = B b (B-spline basis, b ~ N(0, qpriorvar I)),
            % warping coefficients c_i ~ N(0, diag(C)), C_j = vpriorvar / j^vpriordecay,
            % and sigma^2 ~ IG(alpha0, beta0). The warping functions use either
            % the shooting vector (v-transform, |v| < pi) or the centered
            % log-derivative (h-transform, unconstrained) representation.
            %
            % Each sweep updates every c_i by preconditioned (infinite-dimensional)
            % HMC with the exact likelihood gradient, then draws b and sigma^2 from
            % their conjugate full conditionals. Step sizes are adapted per curve
            % during burn-in. The amplitude-phase split is identified only up to a
            % common warp (a ridge that per-curve updates move along slowly), so
            % the centering constraint is imposed by recentering the warps every
            % ncenter iterations so that their Karcher mean is the identity (the
            % template is warped accordingly). The chain is initialized from
            % an elastic (DP) Karcher mean alignment unless init = "identity".
            %
            % Usage:  obj.time_warping_bayes_infHMC()
            %         obj.time_warping_bayes_infHMC(Name=Value)
            %
            % Input (name-value):
            % iter = 2000;          % number of iterations
            % burnin = [];          % number of burn-in iterations (default iter/2)
            % thin = 1;             % thinning of stored posterior draws
            % alpha0 = 0.1;         % inverse gamma prior parameters for sigma^2
            % beta0 = 0.1;
            % h = 0.01;             % initial HMC step size is sqrt(h)
            % L = 10;               % maximum number of leapfrog steps
            % adapt = true;         % adapt step sizes during burn-in
            % target_accept = 0.65; % target acceptance when adapting
            % ncoef = 10;           % number of warping basis coefficients
            % basis = 'fourier';    % 'fourier' or 'legendre' warping basis
            % warp_rep = 'shooting';% 'shooting' (v-transform) or 'clr' (h-transform)
            % vpriorvar = 1;        % prior variance scale for warping coefficients
            % vpriordecay = 1;      % prior variance decays as 1/j^vpriordecay
            % nbasis = 40;          % number of B-spline basis functions for q*
            % qpriorvar = 1e3;      % prior variance of the template coefficients
            % init = "dp";          % "dp" (elastic alignment) or "identity"
            % ncenter = 20;         % recenter every ncenter iterations (0: never)
            % center_all = true;    % recenter after burn-in too (false: burn-in only)
            % extrainfo = true;     % return posterior draws and summaries
            %
            % Output:
            % fdawarp object (fn, qn, q0, fmean, mqn, gam, psi, gamI, stats, qun)
            % with gam the posterior Karcher-type (extrinsic) mean warp of each
            % curve, and obj.mcmc containing
            % mcmc.v_coef: (D,N,nsamp) posterior draws of the warping coefficients
            % mcmc.sigma: sigma draws (all iterations)
            % mcmc.accept: (N,iter-1) HMC acceptance indicators
            % mcmc.epsilon: (N,iter) step sizes
            % mcmc.logl: log-likelihood trace
            % if extrainfo
            % mcmc.gamma_mat: (M,nsamp,N) posterior gammas
            % mcmc.gamma_stats: (M,2,N) pointwise 95% bands of gammas
            % mcmc.q_star_mat: (M,nsamp) posterior templates
            % mcmc.q_star_stats: pointwise 95% bands of the template
            % mcmc.fmean_mat: (M,nsamp) posterior mean functions
            % mcmc.fmean_stats: pointwise 95% bands of the mean function
            arguments
                obj
                mcmcopts.iter = 2000;
                mcmcopts.burnin = [];
                mcmcopts.thin = 1;
                mcmcopts.alpha0 = 0.1;
                mcmcopts.beta0 = 0.1;
                mcmcopts.h = 0.01;
                mcmcopts.L = 10;
                mcmcopts.adapt = true;
                mcmcopts.target_accept = 0.65;
                mcmcopts.ncoef = 10;
                mcmcopts.basis {mustBeMember(mcmcopts.basis,{'fourier','legendre'})} = 'fourier';
                mcmcopts.warp_rep {mustBeMember(mcmcopts.warp_rep,{'shooting','clr'})} = 'shooting';
                mcmcopts.vpriorvar = 1;
                mcmcopts.vpriordecay = 1;
                mcmcopts.nbasis = 40;
                mcmcopts.qpriorvar = 1e3;
                mcmcopts.init {mustBeMember(mcmcopts.init,{'dp','identity'})} = 'dp';
                mcmcopts.ncenter = 20;
                mcmcopts.center_all = true;
                mcmcopts.extrainfo = true;
            end

            [M, N] = size(obj.f);
            iter = mcmcopts.iter;
            burnin = mcmcopts.burnin;
            if isempty(burnin)
                burnin = floor(iter/2);
            end
            if burnin < 1 || burnin >= iter
                error('burnin must be between 1 and iter-1')
            end
            ncoef = mcmcopts.ncoef;
            if strcmpi(mcmcopts.basis,'fourier') && mod(ncoef,2) ~= 0
                error('ncoef must be even for the fourier basis')
            end
            rep = double(strcmpi(mcmcopts.warp_rep,'clr'));

            % work on a uniform grid on [0,1]
            time0 = obj.time(:);
            tn = (time0 - time0(1)) / (time0(end) - time0(1));
            t01 = linspace(0,1,M).';
            f0 = interp1(tn, obj.f, t01, 'linear');
            q = f_to_srvf(f0, t01);

            Phi = bayes_warp_basis(t01, ncoef, mcmcopts.basis);
            C = (mcmcopts.vpriorvar ./ (ceil((1:ncoef)/2).^mcmcopts.vpriordecay)).';
            B = create_basismatrix(t01, mcmcopts.nbasis, 4);
            BtB = B.' * B;
            Bpinv = B \ eye(M);

            % initialization
            if strcmpi(mcmcopts.init,'dp')
                tmp = fdawarp(f0, t01);
                tmp = tmp.time_warping(0, MaxItr=1);
                c_curr = bayes_gam_to_coef(tmp.gam, Phi, rep);
                [gam0, ~] = bayes_coef_to_gam(c_curr, Phi, rep);
                qtil = zeros(M,N);
                for ii = 1:N
                    qtil(:,ii) = warp_q_gamma(q(:,ii), gam0(:,ii), t01);
                end
                qstar = mean(qtil,2);
            else
                c_curr = zeros(ncoef, N);
                dqq = sqrt(sum((q - mean(q,2)).^2,1));
                [~, min_ind] = min(dqq);
                qstar = q(:,min_ind);
            end
            b_curr = Bpinv * qstar;
            qstar = B * b_curr;
            sigma2 = mean(sum((q - qstar).^2,1)) / M;

            % storage
            keep = burnin:mcmcopts.thin:iter;
            nsamp = length(keep);
            v_coef = zeros(ncoef, N, nsamp);
            qstar_mat = zeros(M, nsamp);
            sigma_tr = zeros(1, iter);
            logl_tr = zeros(1, iter);
            accept = zeros(N, iter);
            eps_tr = zeros(N, iter);
            epsilon = sqrt(mcmcopts.h) * ones(N,1);
            da = dual_avg_step(epsilon);
            sigma_tr(1) = sqrt(sigma2);
            eps_tr(:,1) = epsilon;
            ks = 1;

            barobj = ProgressBar(iter-1, 'Title', 'Running MCMC');
            for m = 2:iter
                % update warping coefficients: the curves are conditionally
                % independent given q* and sigma^2, so all are updated in
                % one batched HMC step (one MEX call per leapfrog step)
                [c_curr, aux, accept(:,m), alpha] = group_hmc_update(c_curr, Phi, rep, ...
                    qstar, q, sigma2, C, epsilon, mcmcopts.L);
                qtil = aux(2:end,:);

                % dual-averaging step-size adaptation during burn-in; keep
                % the averaged step sizes afterwards
                if mcmcopts.adapt && m <= burnin
                    [epsilon, da] = dual_avg_step(da, alpha, mcmcopts.target_accept);
                    if m == burnin
                        epsilon = exp(da.log_eps_bar);
                    end
                end

                % template: conjugate Gaussian draw
                P = N*BtB/sigma2 + eye(mcmcopts.nbasis)/mcmcopts.qpriorvar;
                R = chol((P + P.')/2);
                mu = R \ (R.' \ (B.' * sum(qtil,2) / sigma2));
                b_curr = mu + R \ randn(mcmcopts.nbasis,1);
                qstar = B * b_curr;

                % noise variance: conjugate inverse gamma draw
                SSE = sum((qtil - qstar).^2, 1);
                sigma2 = 1/gamrnd(N*M/2 + mcmcopts.alpha0, 1/(sum(SSE)/2 + mcmcopts.beta0));

                % centering constraint (Karcher mean of warps = identity):
                % project onto it every ncenter iterations
                if mcmcopts.ncenter > 0 && (m < burnin || mcmcopts.center_all) && mod(m, mcmcopts.ncenter) == 0
                    gam_c = bayes_coef_to_gam(c_curr, Phi, rep);
                    gam_inv = karcher_mean_inverse(gam_c, t01);
                    for ii = 1:N
                        gam_c(:,ii) = interp1(t01, gam_c(:,ii), gam_inv);
                    end
                    c_curr = bayes_gam_to_coef(gam_c, Phi, rep);
                    qstar = warp_q_gamma(qstar, gam_inv, t01);
                    qstar = qstar(:);
                    b_curr = Bpinv * qstar;
                    qstar = B * b_curr;
                end

                sigma_tr(m) = sqrt(sigma2);
                logl_tr(m) = -sum(SSE)/(2*sigma2) - N*M/2*log(2*pi*sigma2);
                eps_tr(:,m) = epsilon;
                if ks <= nsamp && m == keep(ks)
                    v_coef(:,:,ks) = c_curr;
                    qstar_mat(:,ks) = qstar;
                    ks = ks + 1;
                end
                barobj.step([], [], []);
            end
            barobj.release();

            % posterior summaries of the warping functions
            gamma_mat = zeros(M, nsamp, N);
            gam_mean = zeros(M, N);
            psi_mean = zeros(M, N);
            for ii = 1:N
                [gam_s, psi_s] = bayes_coef_to_gam(squeeze(v_coef(:,ii,:)), Phi, rep);
                gamma_mat(:,:,ii) = interp1(t01, gam_s, tn);
                pm = mean(psi_s, 2);
                pm = pm / sqrt(trapz(t01, pm.^2));
                g = cumtrapz(t01, pm.^2);
                g = (g - g(1)) / (g(end) - g(1));
                gam_mean(:,ii) = interp1(t01, g, tn);
                psi_mean(:,ii) = interp1(t01, pm, tn);
            end

            obj.gam = gam_mean;
            obj.psi = psi_mean;
            obj.fn = zeros(M,N);
            for ii = 1:N
                obj.fn(:,ii) = warp_f_gamma(obj.f(:,ii), obj.gam(:,ii), time0);
            end
            obj.qn = f_to_srvf(obj.fn, time0);
            obj.q0 = f_to_srvf(obj.f, time0);
            % template on the original time scale (SRVF scales by 1/sqrt(T))
            mq01 = mean(qstar_mat, 2);
            obj.mqn = interp1(t01, mq01, tn) / sqrt(time0(end) - time0(1));
            obj.fmean = mean(obj.f(1,:)) + cumtrapz(time0, obj.mqn.*abs(obj.mqn));

            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean, obj.gam(:,ii), time0);
            end
            obj.stats.orig_var = trapz(time0, std_f0.^2);
            obj.stats.amp_var = trapz(time0, std_fn.^2);
            obj.stats.phase_var = trapz(time0, var(fgam,[],2));

            obj.qun = logl_tr(2:end);
            obj.method = "bayesian_infHMC";
            obj.gamI = SqrtMeanInverse(obj.gam);
            obj.rsamps = false;
            obj.type = 'mean';

            obj.mcmc = struct();
            obj.mcmc.v_coef = v_coef;
            obj.mcmc.sigma = sigma_tr;
            obj.mcmc.logl = logl_tr;
            obj.mcmc.accept = accept(:,2:end);
            obj.mcmc.epsilon = eps_tr;
            obj.mcmc.warp_rep = mcmcopts.warp_rep;
            obj.mcmc.burnin = burnin;
            obj.mcmc.iter = iter;
            if mcmcopts.extrainfo
                obj.mcmc.gamma_mat = gamma_mat;
                obj.mcmc.gamma_stats = permute(quantile(gamma_mat, [0.025, 0.975], 2), [1 2 3]);
                obj.mcmc.q_star_mat = interp1(t01, qstar_mat, tn) / sqrt(time0(end) - time0(1));
                obj.mcmc.q_star_stats = quantile(obj.mcmc.q_star_mat, [.025, 0.975], 2);
                obj.mcmc.fmean_mat = srvf_to_f(obj.mcmc.q_star_mat, time0, repmat(mean(obj.f(1,:)), 1, nsamp));
                obj.mcmc.fmean_stats = quantile(obj.mcmc.fmean_mat, [.025, 0.975], 2);
            end
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
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
            % option.penalty = 'roughness'; % warping penalty weighted by lambda (none, roughness, l2gam, l2psi, geodesic)
            % option.spl = true; % use interpolating spline derivative (false: smoothing spline, lossy)
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdawarp object
            arguments
                obj
                lambda = 0;
                option.parallel = 0;
                option.nworkers = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP';
                option.penalty = 'roughness';
                option.spl = true;
                option.MaxItr = 20;
            end
            % time warping on a set of functions
            if option.parallel == 1
                start_parpool(option.nworkers);
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
            q = f_to_srvf(f1,t,~option.spl,option.spl);

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
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                        mf(1), f1(1,k,1), 10, option.penalty);
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq;
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                        mf(1), f1(1,k,1), 10, option.penalty);
                end
            end

            gamI_o = SqrtMeanInverse(gam_o');
            mf = warp_f_gamma(mf,gamI_o,t);
            mq = f_to_srvf(mf,t,~option.spl,option.spl);

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
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                            mf(1,r), f1(1,k,1), 10, option.penalty);
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),t);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t,~option.spl,option.spl);
                        v = q_temp(:,k) - mq_c;
                        d = sqrt(trapz(t, v.*v));
                        if d > eps
                            vtil(:,k) = v/d;
                            dtil(k) = 1/d;
                        end
                    end
                else
                    for k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                            mf(1,r), f1(1,k,1), 10, option.penalty);
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),t);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t,~option.spl,option.spl);
                        v = q_temp(:,k) - mq_c;
                        d = sqrt(trapz(t, v.*v));
                        if d > eps
                            vtil(:,k) = v/d;
                            dtil(k) = 1/d;
                        end
                    end
                end
                q(:,:,r+1) = q_temp;
                f1(:,:,r+1) = f_temp;

                ds_tmp = sqrt(sum(simps(t,(q(:,:,r+1)-mq(:,r)*ones(1,N)).^2)))+lambda*sum(simps(t,(1-sqrt(max(gam_dev',0))).^2));
                if (isreal(ds_tmp))
                    ds(r+1) = ds_tmp;
                else
                    ds(r+1) = abs(ds_tmp);
                end


                % Minimization Step
                % compute the mean of the matched function
                stp = .3;
                if any(dtil)
                    vbar = sum(vtil,2)*sum(dtil)^(-1);
                else
                    vbar = zeros(M,1);
                end
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
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                        mf(1,r), f1(1,k,1), 10, option.penalty);
                end
            else
                for k = 1:N
                    q_c = q(:,k,1); mq_c = mq(:,r);
                    gam_o(k,:) = optimum_reparam(mq_c,q_c,t,lambda,option.method, ...
                        mf(1,r), f1(1,k,1), 10, option.penalty);
                end
            end

            gamI_o = SqrtMeanInverse(gam_o');
            % Center gam and the median template, then warp the original
            % data once by the final gam, so fn = f o gam holds exactly.
            % (Previously fn/qn kept the previous iteration's alignment
            % while gam came from this final matching step.)
            mq(:,r+1) = warp_q_gamma(mq(:,r),gamI_o,obj.time);
            for k = 1:N
                gam_o(k,:) = warp_f_gamma(gam_o(k,:),gamI_o,obj.time);
                f1(:,k,r+1) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                q(:,k,r+1) = f_to_srvf(f1(:,k,r+1),obj.time,~option.spl,option.spl);
            end

            %% Aligned data & stats
            obj.fn = f1(:,:,r+1);
            obj.qn = q(:,:,r+1);
            obj.q0 = q(:,:,1);
            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            obj.mqn = mq(:,r+1);
            obj.fmean = median(obj.f(1,:))+cumtrapz(t,obj.mqn.*abs(obj.mqn));

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
            obj.penalty = option.penalty;
            obj.gamI = gamI_o;
            obj.rsamps = false;
            obj.mcmc = [];
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
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.showplot = 1; % turns on and off plotting
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS, expBayes)
            % option.penalty = 'roughness'; % warping penalty weighted by lambda (none, roughness, l2gam, l2psi, geodesic)
            % option.spl = true; % use interpolating spline derivative (false: smoothing spline, lossy)
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output: structure containing
            % fdawarp object
            arguments
                obj
                mu
                lambda = 0;
                option.parallel = 0;
                option.nworkers = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 1;
                option.method = 'DP';
                option.penalty = 'roughness';
                option.spl = true;
                option.MaxItr = 20;
            end

            % time warping on a set of functions
            if option.parallel == 1
                start_parpool(option.nworkers);
            end

            fprintf('\n lambda = %5.2f \n', lambda);

            [M, N] = size(obj.f);

            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end


            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time,~option.spl,option.spl,option.parallel,option.nworkers);

            %% Compute the q-function of the plot
            mq = f_to_srvf(mu,obj.time,~option.spl,option.spl);

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
                        gam1(k,:) = optimum_reparam(mq,q(:,k),obj.time,lambda,option.method, ...
                            mu(1), obj.f(1,k), 10, option.penalty);
                    end
                    fn1(:,k) = warp_f_gamma(obj.f(:,k,1),gam1(k,:),obj.time);
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time,~option.spl,option.spl);
                end
            else
                for k = 1:N
                    if (strcmpi(option.method,'expBayes'))
                        out_e = pairwise_align_bayes(mu, obj.f(:,k), obj.time, mcmcopts);
                        gam1(k,:) = out_e.gamma;
                    else
                        gam1(k,:) = optimum_reparam(mq,q(:,k),obj.time,lambda,option.method, ...
                            mu(1), obj.f(1,k), 10, option.penalty);
                    end
                    fn1(:,k) = warp_f_gamma(obj.f(:,k,1),gam1(k,:),obj.time);
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time,~option.spl,option.spl);
                end
            end

            obj.gamI = SqrtMeanInverse(gam1');

            %% Aligned data & stats
            obj.q0 = q;
            obj.fn = fn1;
            obj.qn = qn1;
            obj.gam = gam1.';
            std_f0 = std(obj.f, 0, 2);
            std_fn = std(obj.fn, 0, 2);
            obj.mqn = mq;
            obj.fmean = mean(obj.f(1,:))+cumtrapz(obj.time,obj.mqn.*abs(obj.mqn));

            fgam = zeros(M,N);
            for ii = 1:N
                fgam(:,ii) = warp_f_gamma(obj.fmean,obj.gam(:,ii),obj.time);
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
            obj.penalty = option.penalty;
            obj.rsamps = false;
            obj.mcmc = [];
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
            arguments
                obj
                n = 1;
                sort_samples = false;
            end
            if (isempty(obj.type))
                error('Please align first');
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
                ip = zeros(1,n);
                for i = 1:n
                    psi1(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
                    ip(i) = ones(1,M)*psi1(i,:)'/M;
                    len(i) = acos(min(max(ip(i),-1),1));
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
                    len(i) = acos(min(max(ip(i),-1),1));
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

        function [fn, qn, gam] = align_new(obj, f_new, option)
            % ALIGN_NEW align new functions to the fitted template
            % -------------------------------------------------------------------------
            % Aligns new functions to the template srvf (mqn) with the same
            % lambda, method and penalty used to fit the object
            %
            % Usage:  [fn, qn, gam] = obj.align_new(f_new)
            %         [fn, qn, gam] = obj.align_new(f_new, spl=false)
            %
            % Input:
            % f_new: matrix (M,K) of K new functions sampled at obj.time
            % option.spl: srvf derivative used for the fit (default = true);
            %             data smoothed for the fit should be smoothed the
            %             same way before calling
            %
            % Output:
            % fn: (M,K) aligned functions
            % qn: (M,K) aligned srvfs
            % gam: (M,K) warping functions
            arguments
                obj
                f_new
                option.spl = true;
            end
            if (isempty(obj.mqn) || strcmpi(obj.method, 'bayesian'))
                error('Please align fdawarp class with time_warping or time_warping_median first!');
            end
            t = obj.time(:);
            if (isvector(f_new))
                f_new = f_new(:);
            end
            [M, K] = size(f_new);
            if (M ~= length(t))
                error('f_new must have length(obj.time) rows');
            end

            q_new = f_to_srvf(f_new, t, ~option.spl, option.spl);
            fn = zeros(M,K);
            qn = zeros(M,K);
            gam = zeros(M,K);
            for k = 1:K
                gam(:,k) = optimum_reparam(obj.mqn(:), q_new(:,k), t, obj.lambda, ...
                    obj.method, obj.fmean(1), f_new(1,k), 10, obj.penalty);
                fn(:,k) = warp_f_gamma(f_new(:,k), gam(:,k), t);
                qn(:,k) = f_to_srvf(fn(:,k), t, ~option.spl, option.spl);
            end
        end

        function plot(obj)
            % plot plot functional alignment results
            % -------------------------------------------------------------------------
            % Usage: obj.plot()

            % plot random samples only when they are all present,
            % otherwise fall back to the original and aligned data
            has_rsamps = ~isempty(obj.rsamps) && obj.rsamps && ...
                ~isempty(obj.ft) && ~isempty(obj.gams) && ...
                ~isempty(obj.fs) && ~isempty(obj.qs);

            if (has_rsamps)

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
                    plot(obj.time, mean_f0, 'k-', 'linewidth', 1); hold on;
                    plot(obj.time, mean_f0+std_f0, '-', 'linewidth', 1, 'Color', '#66C2A5');
                    plot(obj.time, mean_f0-std_f0, '-', 'linewidth', 1, 'Color', '#66C2A5');
                    title('Original data: Mean \pm STD', 'fontsize', 16);

                    figure(5); clf;
                    plot(obj.time, mean_fn, 'k-', 'linewidth', 1); hold on;
                    plot(obj.time, mean_fn+std_fn, '-', 'linewidth', 1, 'Color', '#66C2A5');
                    plot(obj.time, mean_fn-std_fn, '-', 'linewidth', 1, 'Color', '#66C2A5');
                    title(['Warped data, \lambda = ' num2str(obj.lambda) ': Mean \pm STD'], 'fontsize', 16);

                    figure(6); clf;
                    plot(obj.time, obj.fmean,'LineWidth', 1, 'Color', '#66C2A5');
                    title(['f_{mean}, \lambda = ' num2str(obj.lambda)], 'fontsize', 16);
                end
            end

            % MCMC diagnostics when the Bayesian alignment was used
            if (~isempty(obj.mcmc) && isstruct(obj.mcmc) && ...
                    (isfield(obj.mcmc, 'sigma1') || isfield(obj.mcmc, 'sigma')))
                obj.plot_mcmc();
            end
        end

        function plot_mcmc(obj)
            % PLOT_MCMC plot diagnostics of the Bayesian alignment MCMC chain
            % -------------------------------------------------------------------------
            % Plots trace, posterior and mixing diagnostics of the chain
            % stored in obj.mcmc by time_warping_bayes (pCN) or
            % time_warping_bayes_infHMC. Plots that need
            % the extra information (mcmcopts.extrainfo = true) are only
            % drawn when it is available.
            %
            % Usage: obj.plot_mcmc()
            %
            % Figures:
            %   7 - chain diagnostics: sigma trace, log-likelihood trace,
            %       posterior of sigma, autocorrelation of sigma, running
            %       acceptance rate and pCN step usage (pCN) or HMC
            %       step sizes per function (infHMC)
            %   8 - warping coefficient norm trace for each function
            %   9 - posterior template with 95% credible band
            %  10 - posterior warping functions with 95% credible bands

            if (isempty(obj.mcmc) || ~isstruct(obj.mcmc) || ...
                    ~(isfield(obj.mcmc, 'sigma1') || isfield(obj.mcmc, 'sigma')))
                error('fdawarp:noMCMC', ...
                    'No MCMC output available, run time_warping_bayes or time_warping_bayes_infHMC first');
            end
            mc = obj.mcmc;
            is_hmc = isfield(mc, 'epsilon');

            if (isfield(mc, 'sigma1'))
                sigma1 = mc.sigma1(:);
            else
                sigma1 = mc.sigma(:);
            end
            niter = length(sigma1);
            if (isfield(mc, 'burnin'))
                burnin = mc.burnin;
            elseif (isfield(mc, 'q_star_mat'))
                burnin = niter - size(mc.q_star_mat, 2) + 1;
            else
                burnin = floor(niter/2);
            end
            burnin = min(max(burnin, 1), niter);
            if (isfield(mc, 'logl'))
                logl = mc.logl(:);
            else
                logl = obj.qun(:);
            end
            has_extra = isfield(mc, 'accept') && ~isempty(mc.accept);
            post = burnin:niter;
            col = [0.4 0.76 0.65];

            %% chain diagnostics
            figure(7); clf;
            tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, 'MCMC diagnostics', 'fontsize', 16);

            nexttile;
            plot(1:niter, sigma1, 'Color', col, 'linewidth', 1);
            xline(burnin, 'k--', 'burn-in');
            xlabel('iteration'); ylabel('\sigma');
            title('Trace of \sigma');

            nexttile;
            if (length(logl) == niter)
                plot(1:niter, logl, 'Color', col, 'linewidth', 1);
                xline(burnin, 'k--', 'burn-in');
                xlabel('iteration'); ylabel('log-likelihood');
            else
                axis off;
                text(0.5, 0.5, {'log-likelihood not stored', '(set mcmcopts.extrainfo = true)'}, ...
                    'HorizontalAlignment', 'center');
            end
            title('Trace of log-likelihood');

            nexttile;
            histogram(sigma1(post), 'Normalization', 'pdf', ...
                'FaceColor', col, 'EdgeColor', 'none');
            xline(mean(sigma1(post)), 'k-', 'linewidth', 1.5);
            xline(quantile(sigma1(post), [.025 .975]), 'k--');
            xlabel('\sigma'); ylabel('density');
            title('Posterior of \sigma (95% CI)');

            nexttile;
            maxlag = min(50, length(post)-1);
            rho = acf_vec(sigma1(post), maxlag);
            stem(0:maxlag, rho, 'filled', 'Color', col, 'MarkerSize', 3);
            hold on;
            yline(1.96/sqrt(length(post)) * [-1 1], 'k--');
            hold off;
            xlabel('lag'); ylabel('ACF');
            ess = length(post) / max(1 + 2*sum(rho(2:end)), 1);
            title(sprintf('ACF of \\sigma (ESS \\approx %.0f)', ess));

            nexttile;
            if (has_extra)
                % HMC stores (N,iter-1) indicators: average over functions
                acc = double(mc.accept);
                if (is_hmc)
                    acc = mean(acc, 1);
                    lab = 'warps';
                else
                    lab = 'q^*';
                end
                acc = acc(:);
                it = (1:length(acc)).' + 1;
                plot(it, cumsum(acc)./(1:length(acc)).', 'Color', col, 'linewidth', 1);
                xline(burnin, 'k--', 'burn-in');
                ylim([0 1]);
                xlabel('iteration'); ylabel('rate');
                title(sprintf('Running acceptance of %s (%.2f)', lab, mean(acc)));
            else
                axis off;
                text(0.5, 0.5, {'acceptance not stored', '(set mcmcopts.extrainfo = true)'}, ...
                    'HorizontalAlignment', 'center');
                title('Running acceptance of q^*');
            end

            nexttile;
            if (is_hmc)
                eps_tr = mc.epsilon;
                semilogy(1:size(eps_tr, 2), eps_tr.', 'linewidth', 1);
                xline(burnin, 'k--', 'burn-in');
                xlabel('iteration'); ylabel('\epsilon');
                title('HMC step size per function');
            elseif (has_extra && isfield(mc, 'betas_ind'))
                bind = mc.betas_ind(:);
                acc = logical(mc.accept(:));
                if (isfield(mc, 'zpcn'))
                    nb = length(mc.zpcn.betas);
                else
                    nb = max(bind);
                end
                nprop = accumarray(bind(bind>0), 1, [nb 1]);
                nacc = accumarray(bind(bind>0 & acc), 1, [nb 1]);
                bar(1:nb, [nprop - nacc, nacc], 'stacked', 'EdgeColor', 'none');
                if (isfield(mc, 'zpcn'))
                    xticklabels(compose('%g', mc.zpcn.betas));
                    xlabel('pCN \beta');
                else
                    xlabel('pCN step index');
                end
                ylabel('proposals');
                legend({'rejected', 'accepted'}, 'Location', 'best');
                rate = nacc ./ max(nprop, 1);
                text(1:nb, nprop, compose('%.2f', rate), 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom');
                ylim([0 1.15*max([nprop; 1])]);
                title('q^* proposals by pCN step (acc. rate)');
            else
                axis off;
                text(0.5, 0.5, {'pCN step usage not stored', '(set mcmcopts.extrainfo = true)'}, ...
                    'HorizontalAlignment', 'center');
                title('pCN step usage');
            end

            %% warping coefficient traces
            if (isfield(mc, 'g_coef') && ~isempty(mc.g_coef))
                figure(8); clf;
                gnorm = squeeze(sqrt(sum(mc.g_coef.^2, 2)));
                if (isvector(gnorm))
                    gnorm = gnorm(:);
                end
                plot(1:size(gnorm, 1), gnorm, 'linewidth', 1);
                xline(burnin, 'k--', 'burn-in');
                xlabel('iteration'); ylabel('norm of warping coefficients');
                title('Trace of warping coefficient norm per function', 'fontsize', 16);
            elseif (isfield(mc, 'v_coef') && ~isempty(mc.v_coef))
                % infHMC stores (D,N,nsamp) post burn-in draws only
                figure(8); clf;
                vnorm = reshape(sqrt(sum(mc.v_coef.^2, 1)), size(mc.v_coef, 2), []).';
                plot(1:size(vnorm, 1), vnorm, 'linewidth', 1);
                xlabel('posterior draw (after burn-in)'); ylabel('norm of warping coefficients');
                title('Trace of warping coefficient norm per function', 'fontsize', 16);
            end

            if (~isfield(mc, 'fmean_mat') || ~isfield(mc, 'gamma_mat'))
                return
            end

            %% posterior template
            figure(9); clf;
            ts = linspace(obj.time(1), obj.time(end), size(mc.fmean_mat, 1)).';
            ndraw = size(mc.fmean_mat, 2);
            idx = unique(round(linspace(1, ndraw, min(ndraw, 100))));
            plot(ts, mc.fmean_mat(:, idx), 'Color', [0.7 0.7 0.7 0.3], 'linewidth', 0.5);
            hold on;
            fill([ts; flipud(ts)], [mc.fmean_stats(:,1); flipud(mc.fmean_stats(:,2))], ...
                col, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            h = plot(obj.time, obj.fmean, 'k-', 'linewidth', 2);
            hold off;
            legend(h, 'posterior mean', 'Location', 'best');
            title('Posterior template with 95% credible band', 'fontsize', 16);

            %% posterior warping functions
            figure(10); clf;
            N = size(mc.gamma_mat, 3);
            tg = linspace(0, 1, size(mc.gamma_mat, 1)).';
            colors = lines(N);
            hold on;
            for jj = 1:N
                band = quantile(mc.gamma_mat(:,:,jj), [.025 .975], 2);
                fill([tg; flipud(tg)], [band(:,1); flipud(band(:,2))], colors(jj,:), ...
                    'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            for jj = 1:N
                plot(tg, obj.gam(:,jj), 'Color', colors(jj,:), 'linewidth', 1);
            end
            plot([0 1], [0 1], 'k:');
            hold off;
            axis square; box on;
            title('Posterior warping functions with 95% credible bands', 'fontsize', 16);
        end

    end

end

function rho = acf_vec(x, maxlag)
% sample autocorrelation of x for lags 0:maxlag
x = x(:) - mean(x);
n = length(x);
c0 = sum(x.^2);
rho = zeros(maxlag+1, 1);
if (c0 == 0)
    rho(1) = 1;
    return
end
for k = 0:maxlag
    rho(k+1) = sum(x(1:n-k).*x(1+k:n)) / c0;
end
end

function out = statsFun(vec)
a = quantile(vec,0.025);
b = quantile(vec,0.975);
out = [a,b];
end

function out = f_exp1(g)
out.x = g.x;
out.y = bcalcY(f_L2norm(g), g.y);
end

function [x,yy] = f_exp1inv(psi)
x = psi.x;
y = psi.y;

[x,ia,~] = unique(x);
[x, ia1] = sort(x);
y = y(ia);
y = y(ia1);
inner = round(trapzCpp(x,y), 10);

if ((inner < 1.001) && (inner >= 1))
    inner = 1;
end
if ((inner <= -1) && (inner > -1.001))
    inner = -1;
end
if ((inner < (-1)) || (inner > 1))
    fprintf("exp1inv: can't calculate the acos of: %f\n", inner);
end

theta = acos(inner);
yy = theta / sin(theta) .* (y - repelem(inner,length(y))');

if (theta==0)
    yy = zeros(1,length(x));
end

end

% function for calculating the next MCMC sample given current state
function [g_coef, accept, zpcnInd] = f_updateg_pw(g_coef_curr,g_basis,var1_curr,q1,q2,propose_g_coef)
g_coef_prop = propose_g_coef(g_coef_curr);

tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
while (min(tst.y)<0)
    g_coef_prop = propose_g_coef(g_coef_curr);
    tst = f_exp1(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis, false));
end

SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis, false), q1, q2);

SSE_prop = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_prop.prop,g_basis,false), q1, q2);

logl_ratio = -(SSE_prop - SSE_curr)/(2*var1_curr);

ratio = min(1, exp(logl_ratio));

u = rand;
if (u <= ratio)
    g_coef = g_coef_prop.prop;
    accept = true;
    zpcnInd = g_coef_prop.ind;
end

if (u > ratio)
    g_coef = g_coef_curr;
    accept = false;
    zpcnInd = g_coef_prop.ind;
end
end

% function for calculating the next MCMC sample given current state
function [q_star, accept, zpcnInd] = f_updateq(g_coef_curr,g_basis,var1_curr,q,q_star_coef_curr,q_basis,SSE_curr,propose_q_star)
q_star_coef_prop = propose_q_star(q_star_coef_curr);
ind = q_star_coef_prop.ind;
q_star_prop = f_basistofunction(q_basis.x,0,q_star_coef_prop.prop,q_basis, false);
q_star_curr = f_basistofunction(q_basis.x,0,q_star_coef_curr,q_basis, false);

if (SSE_curr == 0)
    SSE_curr = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis, false), q, q_star_curr);
end

SSE_prop = f_SSEg_pw(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_prop.y);

logl_curr = f_logl(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_curr.y, var1_curr, SSE_curr);

logl_prop = f_logl(f_basistofunction(g_basis.x,0,g_coef_curr,g_basis,false), q, q_star_prop.y, var1_curr, SSE_prop);

ratio = min(1, exp(logl_prop-logl_curr));

u = rand;
if (u <= ratio)
    q_star = q_star_coef_prop.prop;
    accept = true;
    zpcnInd = ind;
end

if (u > ratio)
    q_star = q_star_coef_curr;
    accept = false;
    zpcnInd = ind;
end
end

%##########################################################################
% For pairwise registration, evaluate the loglikelihood of g, given q1 and
% q2
% g, q1, q2 are all given in the form of struct.x and struct.y
% var1: model variance
% SSEg: if not provided, than re-calculate
% returns a numeric value which is logl(g|q1,q2), see Eq 10 of JCGS
%##########################################################################
% SSEg: sum of sq error= sum over ti of { q1(ti)-{q2,g}(ti) }^2
% (Eq 11 of JCGS)
function out = f_SSEg_pw(g, q, q_star)
obs_domain = g.x;
out = zeros(1,size(q.y,2));
for ii = 1:size(q.y,2)
    g_tmp = g;
    g_tmp.y = g_tmp.y(:,ii);
    exp1g_temp = f_predictfunction(f_exp1(g_tmp), obs_domain, 0);
    pt = [0; bcuL2norm2(obs_domain, exp1g_temp.y)];
    q_tmp = q;
    q_tmp.y = q_tmp.y(:,ii);
    tmp = f_predictfunction(q_tmp, pt, 0);
    vec = (tmp.y .* exp1g_temp.y - q_star).^2;
    out(ii) = sum(vec);
end
end

function out = f_logl(g, q, q_star, var1, SSEg)
if (SSEg == 0)
    SSEg = f_SSEg_pw(g, q, q_star);
end
n = length(q.x);
out = n * log(1/sqrt(2*pi)) - n * log(sqrt(var1)) - (SSEg ./ (2 * var1));
out = sum(out);
end

%##########################################################################
% Extrapolate a function given by a discrete vector
% f: function in the form of list$x, list$y
% at: t values for which f(t) is returned
% deriv: can calculate derivative
% method: smoothing method: 'linear' (default), 'cubic'
% returns: $y==f(at), $x==at
%##########################################################################
function out = f_predictfunction(f, at, deriv)
if (deriv == 0)
    out.y = zeros(length(at),size(f.y,2));
    out.x = at;
    for ii = 1:size(f.y,2)
        result = interp1(f.x,f.y(:,ii),at,'linear','extrap');
        out.y(:,ii) = result;
    end
end

if (deriv == 1)
    fmod = interp1(f.x,f.y,at,'linear','extrap');
    diffy1 = [0; diff(fmod)];
    diffy2 = [diff(fmod); 0];
    diffx1 = [0; diff(at)];
    diffx2 = [diff(at); 0];


    out.x = at;
    out.y = (diffy2 + diffy1) ./ (diffx2 + diffx1);
end
end

%##########################################################################
% calculate L2 norm of a function, using trapezoid rule for integration
% f:function in the form of list$x, list$y
% returns ||f||, a numeric value
%##########################################################################
function out = f_L2norm(f)
out = border_l2norm(f.x,f.y);
end

%##########################################################################
% Different basis functions b_i()
% f.domain: grid on which b_i() is to be returned
% numBasis: numeric value, number of basis functions used
% (note: #basis = #coef/2 for Fourier basis)
% fourier.p: period of the Fourier basis used
% returns a struct:
%     matrix: with nrow=length(t) and ncol=numBasis (or numBasis*2 for
%     Fourier)
%     x: f.domain
%##########################################################################
function out = basis_fourier(f_domain, numBasis, fourier_p)
result = zeros(length(f_domain), 2*numBasis);
for i = 1:(2*numBasis)
    j = ceil(i/2);
    if (mod(i,2) == 1)
        result(:,i) = sqrt(2) * sin(2*j*pi*f_domain./fourier_p);
    end
    if (mod(i,2) == 0)
        result(:,i) = sqrt(2) * cos(2*j*pi*f_domain./fourier_p);
    end
    out.x = f_domain;
    out.matrix = result;
end
end

function out = basis_bspline(f_domain, numBasis, df)
out.x = f_domain;
out.matrix = create_basismatrix(f_domain,numBasis,df);
end

%##########################################################################
% Given the coefficients of basis functions, returns the actual function on
% a grid
% f.domain: numeric vector, grid of the actual function to return
% coefconst: leading constant term
% coef: numeric vector, coefficients of the basis functions
%       Note: if #coef < #basis functions, only the first %coef basis
%             functions will be used
% basis: in the form of list$x, list$matrix
% plotf: if true, show a plot of the function generated
% returns the generated function in the form of struct.x=f.domain, struct.y
%##########################################################################
function result = f_basistofunction(f_domain, coefconst, coef, basis, plotf)
if (size(basis.matrix,2) < size(coef,1))
    error('In f_basistofunction, #coeffients exceeds #basis functions.')
end
result.x = basis.x;
result.y = basis.matrix(:,(1:size(coef,1))) * coef + coefconst;

result = f_predictfunction(result, f_domain, 0);
if (plotf)
    plot(result.x,result.y)
end
end

%##########################################################################
% Calculate Karcher mean/median with Alg1/Alg2 in (Kurtek,2014)
% x: vector of length = length(domain of the psi's)
% y: M columns, each of length = length(x)
% e1, e2: small positive constants
% method: 'ext' = extrinsic, 'int' = intrinsic
% returns posterier mean/median of M psi's (of form .x, .y)
%##########################################################################
function out = f_psimean(x, y)
rmy = mean(y,2);
tmp.x = x;
tmp.y = rmy;
result = rmy / f_L2norm(tmp);
out.x = x;
out.y = result;
end

%##########################################################################
% calculate phi(gamma), phiinv(psi)
% gamma, psi: function in the form of struct.x, struct.y
% returns phi(gamma) or phiinv(psi), function in the form of struct.x,
% struct.y
%##########################################################################
function result = f_phi(gamma)
f_domain = gamma.x;
k = f_predictfunction(gamma, f_domain, 1);
k = k.y;
if (isempty(find(k < 0, 1)) ~= 0)
    idx = k < 0;
    k(idx) = 0;
end
result.x = f_domain;
result.y = sqrt(k);
if (f_L2norm(result) >= (1.01) || f_L2norm(result) <= (0.99))
    result.y = result.y / f_L2norm(result);
end
end
% the function returned by phi(gamma) = psi is always positive and has L2norm 1

function out = f_phiinv(psi)
f_domain = psi.x;
result = [0; bcuL2norm2(f_domain, psi.y)];
out.x = f_domain;
out.y = result;
end

% Normalize gamma to [0,1]
function gam = norm_gam(gam)
gam = (gam-gam(1))./(gam(end)-gam(1));
end

% geometry of the group alignment likelihood for a batch of curves (one per
% column): aux = [SSE; q_i warped]
function [nll, g, aux] = group_geom(c, Phi, rep, qstar, Q, sigma2)
[nll, g, SSE, ~, q2w] = bayes_warp_geom(c, Phi, rep, qstar, Q, sigma2);
aux = [SSE; q2w];
end

% one batched infinite-dimensional HMC update of the warps of curves Q
function [c, aux, accept, alpha] = group_hmc_update(c, Phi, rep, qstar, Q, sigma2, C, epsilon, L)
geom = @(x) group_geom(x, Phi, rep, qstar, Q, sigma2);
[nll, g, aux] = geom(c);
[c, ~, ~, aux, accept, alpha] = infhmc_step_batch(c, nll, g, aux, geom, C, epsilon, L);
end

% Inverse of the Karcher mean of warping functions gam (T,N) on the uniform
% grid t; same estimate as SqrtMeanInverse, vectorized over the curves and
% started from the extrinsic mean (the warps are near identity after the
% first centering step, so a few full gradient steps converge)
function gamI = karcher_mean_inverse(gam, t)
[T, N] = size(gam);
dt = t(2) - t(1);
w = dt * ones(T, 1);
w([1, T]) = dt / 2;
[~, gd] = gradient(gam, 1, dt); % derivative along each column
psi = sqrt(max(gd, 0));
psi = psi ./ sqrt(w.' * psi.^2);
mu = mean(psi, 2);
mu = mu / sqrt(w.' * mu.^2);
for it = 1:100
    ip = min(max(w.' * (psi .* mu), -1), 1);
    th = acos(ip);
    fac = ones(1, N);
    nz = th > 1e-10;
    fac(nz) = th(nz) ./ sin(th(nz));
    vbar = mean((psi - mu .* ip) .* fac, 2);
    nv = sqrt(w.' * vbar.^2);
    if nv < 1e-8
        break
    end
    mu = cos(nv) * mu + sin(nv) / nv * vbar;
    mu = mu / sqrt(w.' * mu.^2);
end
gam_mu = cumtrapz(t, mu.^2);
gam_mu = (gam_mu - gam_mu(1)) / (gam_mu(end) - gam_mu(1));
gamI = invertGamma(gam_mu);
gamI = gamI(:);
end
