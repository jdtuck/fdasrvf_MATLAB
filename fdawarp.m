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
    %   time_warping_bayes - align functions and find karcher mean using Bayesian
    %   approach
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

            if (size(f,1) ~= length(time))
                error('Columns of f and time must be equal')
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
            % option.spl = true; % use spline interpolation 
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
                option.spl = true;
                option.MaxItr = 20;
            elseif nargin < 3
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
                option.spl = true;
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
            q = f_to_srvf(f1,obj.time,option.spl);

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
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time,option.spl);
                    end
                else
                    for k = 1:N
                        q_c = q(:,k,1); mq_c = mq(:,r);
                        gam_o(k,:) = optimum_reparam(mq_c,q_c,obj.time,lambda,option.method,option.w, ...
                            mf(1,r), f1(1,k,1));
                        gam_dev(k,:) = gradient(gam_o(k,:), 1/(M-1));
                        f_temp(:,k) = warp_f_gamma(f1(:,k,1),gam_o(k,:),obj.time);
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time,option.spl);
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
            obj.gamI = gamI_o;
            obj.rsamps = false;
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
            [M, N] = size(obj.f);

            if nargin < 2
                mcmcopts.iter = 2000;
                mcmcopts.burnin = min(5e3,mcmcopts.iter/2);
                mcmcopts.alpha0 = 0.1;
                mcmcopts.beta0 = 0.1;
                tmp.betas = [0.5,0.5,0.005,0.0001];
                tmp.probs = [0.1,0.1,0.7,0.1];
                mcmcopts.zpcn = tmp;
                mcmcopts.propvar = 4;
                mcmcopts.initcoef = zeros(20, N);
                mcmcopts.nbasis = 30;
                mcmcopts.npoints = 200;
                mcmcopts.extrainfo = true;
                mcmcopts.ncenter = 100;
            end

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
%             q_basis.matrix = [ones(numSimPoints,1) q_basis.matrix];
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
                sdvec = [pw_sim_global_sigma_q./repelem(1:pw_sim_global_Mq,2)];
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
                    gamI_o = SqrtMeanInverse(result_posterior_gamma.');
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

            obj.mcmc.g_coef = result.g_coef;
            obj.mcmc.sigma1 = result.sigma1;

            if (mcmcopts.extrainfo)
                obj.mcmc.accept = result.accept(2:end);
                obj.mcmc.betas_ind = result.accept_betas(2:end);
                obj.mcmc.gamma_mat = gamma_mat;
                obj.mcmc.gamma_stats = gamma_stats;
                obj.mcmc.q_star_mat = result.q_star(valid_index,:)';
                obj.mcmc.q_star_stats = quantile(obj.mcmc.q_star_mat,[.025,0.975],2);
                obj.mcmc.fmean_mat = srvf_to_f(obj.mcmc.q_star_mat,q_basis.x,repelem(mean(obj.f(1,:)),length(valid_index)));
                obj.mcmc.fmean_stats = quantile(obj.mcmc.fmean_mat,[.025,0.975],2);
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
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, DP2, SIMUL, RBFGS)
            % option.w = 0.0; % BFGS weight
            % option.spl = true; % use spline interpolation
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
                option.spl = true;
                option.MaxItr = 20;
            elseif nargin < 3
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.w = 0.0;
                option.spl = true;
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
            q = f_to_srvf(obj.f,t,option.spl);

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
            mq = f_to_srvf(mf,t,option.spl);

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
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t,option.spl);
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
                        q_temp(:,k) = f_to_srvf(f_temp(:,k),t,option.spl);
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
            % option.spl = true; % use spline interpolation
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
                option.spl = true;
                option.MaxItr = 20;
            elseif nargin < 4
                option.parallel = 0;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.showplot = 1;
                option.method = 'DP1';
                option.w = 0.0;
                option.spl = true;
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
            q = f_to_srvf(obj.f,obj.time,option.spl);

            %% Compute the q-function of the plot
            mq = f_to_srvf(mu,obj.time,option.spl);

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
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time,option.spl);
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
                    qn1(:,k) = f_to_srvf(fn1(:,k),obj.time,option.spl);
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
