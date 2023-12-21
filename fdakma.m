classdef fdakma
    %fdakma A class to provide a kmeans clustering and alignment
    % -------------------------------------------------------------------------
    % This class provides kmeans and alignment using the
    % SRVF framework
    %
    % Usage:  obj = fdakma(f,t)
    %
    % where:
    %   f: (M,N): matrix defining N functions of M samples
    %   time: time vector of length M
    %
    %
    % fdakma Properties:
    %   f - original functions
    %   time - time
    %   fn - aligned functions in cell of clusters
    %   qn - aligned srvfs in cell of clusters
    %   gam - warping functions in cell of clusters
    %   gamI - inverse gamma
    %   q0 - original srvfs
    %   labels - cluster labels
    %   templates - cluster centers
    %   templates_q - cluster srvf centers
    %   qun - cost function
    %   method - optimization method
    %
    %
    % fdakma Methods:
    %   fdakma - class constructor
    %   kmeans - perform kmeans
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        f            % original functions
        time         % time
        fn           % aligned functions in cell of clusters
        qn           % aligned srvfs in cell of clusters
        gam          % warping functions in cell of clusters
        gamI         % inverse gamma
        q0           % original srvfs
        labels       % cluster labels
        templates    % cluster centers
        templates_q  % cluster srvf centers
        qun          % cost function
        lambda       % lambda
        method       % optimization method
    end
    
    methods
        function obj = fdakma(f,time)
            %fdakma Construct an instance of this class
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
        
        function obj = kmeans(obj,K,seeds,lambda1, option)
            % KMEANS K-Means clustering and alignment
            % -------------------------------------------------------------------------
            % This function clusters functions and aligns using the elastic square-root
            % slope (srsf) framework.
            %
            % Usage:  obj.kmeans(K)
            %         obj.kmeans(K,seeds)
            %         obj.kmeans(K,seeds,lambda1)
            %
            % Input:
            % K: number of clusters
            % seeds: indexes of cluster center functions (default [])
            % lambda1: controls amount of alignment
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.alignment = true; % performa alignment
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
            % option.MaxItr = 20;  % maximum iterations
            % option.thresh = 0.01; % cost function threshold
            %
            % Output:
            % fdakma object
            
            if nargin < 2
                seeds = [];
                obj.lambda = 0;
                option.parallel = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.MaxItr = 20;
                option.thresh = 0.01;
            elseif nargin < 3
                obj.lambda = 0;
                option.parallel = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.MaxItr = 20;
                option.thresh = 0.01;
            elseif nargin < 4
                obj.lambda = lambda1;
                option.parallel = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.MaxItr = 20;
                option.thresh = 0.01;
            else
                obj.lambda = lambda1;
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
            
            fprintf('\n lambda = %5.1f \n', obj.lambda);
            
            
            [M, N] = size(obj.f);
            
            if (isempty(seeds))
                template_ind = randsample(N,K);
            else
                template_ind = seeds;
            end
            obj.templates = zeros(M,K);
            for i=1:K
                obj.templates(:,i) = obj.f(:,template_ind(i));
            end
            cluster_id = zeros(1,N);
            obj.qun = zeros(0, option.MaxItr);
            
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            
            
            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time);
            obj.q0 = q;
            obj.templates_q = zeros(M,K);
            for i=1:K
                obj.templates_q(:,i) = q(:,template_ind(i));
            end
            
            %% Set initial using the original f space
            fprintf('\nInitializing...\n');
            f_temp = zeros(length(obj.time),N);
            q_temp = zeros(length(obj.time),N);
            for r = 1:option.MaxItr
                gam1 = {};
                Dy = zeros(K,N);
                qn1 = {};
                fn1 = {};
                fprintf('updating step: r=%d\n', r);
                if r == option.MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end
                for i=1:K
                    gamt = zeros(N,size(q,1));
                    if option.parallel == 1
                        parfor k = 1:N
                            q_c = q(:,k); mq_c = obj.templates_q(:,i);
                            if (option.alignment)
                                gamt(k,:) = optimum_reparam(mq_c,q_c,obj.time,obj.lambda,option.method, ...
                                    obj.templates(1,i), obj.f(1,k));
                            else
                                gamt(k,:) = linspace(0,1,M);
                            end
                            f_temp(:,k) = warp_f_gamma(obj.f(:,k),gamt(k,:),obj.time);
                            q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                            Dy(i,k) = sqrt(sum(trapz(obj.time,(q_temp(:,k)-obj.templates_q(:,i)).^2)));
                        end
                    else
                        for k = 1:N
                            q_c = q(:,k); mq_c = obj.templates_q(:,i);
                            if (option.alignment)
                                gamt(k,:) = optimum_reparam(mq_c,q_c,obj.time,obj.lambda,option.method, ...
                                    obj.templates(1,i), obj.f(1,k));
                            else
                                gamt(k,:) = linspace(0,1,M);
                            end
                            f_temp(:,k) = warp_f_gamma(obj.f(:,k),gamt(k,:),obj.time);
                            q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                            Dy(i,k) = sqrt(sum(trapz(obj.time,(q_temp(:,k)-obj.templates_q(:,i)).^2)));
                        end
                    end
                    gam1{i} = gamt;
                    qn1{i} = q_temp;
                    fn1{i} = f_temp;
                end
                [~,cluster_id] = min(Dy,[],1);
                
                %% Normalization
                for i=1:K
                    id = cluster_id == i;
                    ftmp = fn1{i}(:,id);
                    gamtmp = gam1{i}(id,:);
                    obj.gamI = SqrtMeanInverse(gamtmp);
                    N1 = size(ftmp,2);
                    f_temp1 = zeros(M,N1);
                    q_temp1 = zeros(M,N1);
                    gam2 = zeros(N1,M);
                    if option.parallel == 1
                        parfor k = 1:N1
                            f_temp1(:,k) = warp_f_gamma(ftmp(:,k),obj.gamI,obj.time);
                            q_temp1(:,k) = f_to_srvf(f_temp1(:,k),obj.time);
                            gam2(k,:) = warp_f_gamma(gamtmp(k,:),obj.gamI,obj.time);
                        end
                    else
                        for k = 1:N1
                            f_temp1(:,k) = warp_f_gamma(ftmp(:,k),obj.gamI,obj.time);
                            q_temp1(:,k) = f_to_srvf(f_temp1(:,k),obj.time);
                            gam2(k,:) = warp_f_gamma(gamtmp(k,:),obj.gamI,obj.time);
                        end
                    end
                    gam1{i}(id,:) = gam2;
                    qn1{i}(:,id) = q_temp1;
                    fn1{i}(:,id) = f_temp1;
                end
                
                %% Template Identification
                qun_t = zeros(1,K);
                for i = 1:K
                    id = cluster_id == i;
                    old_templates_q = obj.templates_q;
                    obj.templates_q(:,i) = mean(qn1{i}(:,id),2);
                    obj.templates(:,i) = mean(fn1{i}(:,id),2);
                    qun_t(i) = norm(obj.templates_q(:,i)-old_templates_q(:,i))/norm(old_templates_q(:,i));
                end
                obj.qun(r) = mean(qun_t);
                
                if obj.qun(r) < option.thresh || r >= option.MaxItr
                    break;
                end
            end
            
            %% Output
            ftmp = {};
            qtmp = {};
            gamtmp = {};
            for i = 1:K
                id = cluster_id == i;
                ftmp{i} = fn1{i}(:,id);
                qtmp{i} = qn1{i}(:,id);
                gamtmp{i} = gam1{i}(id,:);
            end
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            obj.fn = ftmp;
            obj.qn = qtmp;
            obj.labels = cluster_id;
            obj.method = option.method;
            obj.qun = obj.qun(1:r);
        end
        
        function obj = DPmeans(obj,lambda1,option)
            % DPmeans Dirichlet Process K-Means clustering and alignment
            % -------------------------------------------------------------------------
            % This function clusters functions and aligns using the elastic square-root
            % slope (srsf) framework using Dirchelet Process Mixture Model
            %
            % Usage:  obj.kmeans(lambda)
            %         obj.kmeans(lambda,option)
            %
            % Input:
            % lambda1 : cluster penalty parameter
            %
            % default options
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.alignment = true; % performa alignment
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP1'; % optimization method (DP, SIMUL, RBFGS)
            % option.MaxItr = 20;  % maximum iterations
            % option.thresh = 0.01; % cost function threshold
            %
            % Output:
            % fdakma object
            
            if nargin < 2
                obj.lambda = NaN;
                option.parallel = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.MaxItr = 20;
                option.thresh = 0.01;
            elseif nargin < 3
                obj.lambda = lambda1;
                option.parallel = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP1';
                option.MaxItr = 20;
                option.thresh = 0.01;
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
            [M, N] = size(obj.f);
            cluster_id = ones(1,N);
            obj.qun = zeros(0, option.MaxItr);
            
            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time);
            obj.q0 = q;
            
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            mnq = mean(q,2);
            dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
            [~, min_ind] = min(dqq);
            K = 1;
            obj.templates(:,1) = obj.f(:,min_ind);
            obj.templates_q(:,1) = q(:,min_ind);

            if (isnan(obj.lambda))
                fprintf('\n Initializing Lambda \n');
                Kinit = 4;
                obj.lambda = kpp_init(obj.f,q,obj.time,Kinit,0,option);
            end
            
            fprintf('\n lambda = %5.1f \n', obj.lambda);
            
            %% Set initial using the original f space
            fprintf('\nInitializing...\n');
            f_temp = zeros(length(obj.time),N);
            q_temp = zeros(length(obj.time),N);
            for r = 1:option.MaxItr
                gam1 = {};
                Dy = zeros(K,N);
                qn1 = {};
                fn1 = {};
                fprintf('updating step: r=%d\n', r);
                if r == option.MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end
                for i=1:K
                    gamt = zeros(N,size(q,1));
                    if option.parallel == 1
                        parfor k = 1:N
                            q_c = q(:,k); mq_c = obj.templates_q(:,i);
                            if (option.alignment)
                                gamt(k,:) = optimum_reparam(mq_c,q_c,obj.time,obj.lambda,option.method, ...
                                    obj.templates(1,i), obj.f(1,k));
                            else
                                gamt(k,:) = linspace(0,1,M);
                            end
                            f_temp(:,k) = warp_f_gamma(obj.f(:,k),gamt(k,:),obj.time);
                            q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                            Dy(i,k) = sqrt(sum(trapz(obj.time,(q_temp(:,k)-obj.templates_q(:,i)).^2)));
                        end
                    else
                        for k = 1:N
                            q_c = q(:,k); mq_c = obj.templates_q(:,i);
                            if (option.alignment)
                                gamt(k,:) = optimum_reparam(mq_c,q_c,obj.time,obj.lambda,option.method, ...
                                    obj.templates(1,i), obj.f(1,k));
                            else
                                gamt(k,:) = linspace(0,1,M);
                            end
                            f_temp(:,k) = warp_f_gamma(obj.f(:,k),gamt(k,:),obj.time);
                            q_temp(:,k) = f_to_srvf(f_temp(:,k),obj.time);
                            Dy(i,k) = sqrt(sum(trapz(obj.time,(q_temp(:,k)-obj.templates_q(:,i)).^2)));
                        end
                    end
                    gam1{i} = gamt;
                    qn1{i} = q_temp;
                    fn1{i} = f_temp;
                end
                
                %% Template Identification
                [minval,cluster_id] = min(Dy,[],1);
                idx = minval > obj.lambda;
                for ii = 1:length(idx)
                    if idx(ii) ~= 0
                        K = K + 1;
                        cluster_id(ii) = K;
                        obj.templates(:,K) = obj.f(:,ii);
                        obj.templates_q(:,K) = q(:,ii);
                    end
                end
                
                qun1 = zeros(1,K);
                for i = 1:K
                    id = cluster_id == i;
                    for j = 1:length(id)
                        qun1(i) = qun1(i) + sqrt(sum(trapz(obj.time,(q(:,id(j))-obj.templates_q(:,i)).^2)));
                    end
                end
                obj.qun(r+1) = sum(qun1)+obj.lambda*K;
                
                if abs(obj.qun(r+1)-obj.qun(r)) < option.thresh || r >= option.MaxItr
                    break;
                end
            end
            
            %% Output
            ftmp = {};
            qtmp = {};
            gamtmp = {};
            for i = 1:K
                id = cluster_id == i;
                ftmp{i} = fn1{i}(:,id);
                qtmp{i} = qn1{i}(:,id);
                gamtmp{i} = gam1{i}(id,:);
            end
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            obj.fn = ftmp;
            obj.qn = qtmp;
            obj.labels = cluster_id;
            obj.method = option.method;
            obj.qun = obj.qun(1:r);
        end
        
        function plot(obj)
            % PLOT plot kmeans clustering alignment results
            % -------------------------------------------------------------
            % Usage: obj.plot()
            figure(1); clf;
            plot(obj.time, obj.f, 'linewidth', 1);
            title('Original data', 'fontsize', 16);
            
            K = length(obj.fn);
            colors = varycolor(K);
            figure(2); clf; hold all
            for k=1:K
                plot(obj.time, obj.templates(:,k), 'Color', colors(k,:));
            end
            title('Cluster Mean Functions');
            
            figure(3); clf;
            plot(obj.time, obj.fn{1}, 'Color', colors(1,:));
            hold all
            for k = 2:K
                plot(obj.time, obj.fn{k}, 'Color', colors(k,:));
            end
            for k = 1:K
                plot(obj.time, obj.templates(:,k), 'Color', colors(k,:));
            end
            title('Clustered Functions');
            
            figure(4); clf;
            plot(obj.time, obj.qn{1}, 'Color', colors(1,:));
            hold all
            for k = 2:K
                plot(obj.time, obj.qn{k}, 'Color', colors(k,:));
            end
            for k = 1:K
                plot(obj.time, obj.templates_q(:,k), 'Color', colors(k,:));
            end
            title('Clustered Functions');
        end
    end
end

function [lambda] = kpp_init(f,q,time,K,lambda,option)
%k++ init
%lambda: max distance to k++ means
[M,N] = size(q);
mu = zeros(M,K);
mu_f = mu;
idx = ceil(rand*N);
mu(:,1) = q(:,idx);
mu_f(:,1) = f(:,idx);
Dy = inf(K,N);
for i = 2:K
    if option.parallel == 1
        parfor k = 1:N
            q_c = q(:,k); mq_c = mu(:,i-1);
            if (option.alignment)
                gamt = optimum_reparam(mq_c,q_c,time,lambda,option.method, ...
                    mu_f(1,i-1), f(1,k));
            else
                gamt = linspace(0,1,M);
            end
            f_temp = warp_f_gamma(f(:,k),gamt,time);
            q_temp = f_to_srvf(f_temp,time);
            Dy(i,k) = sqrt(sum(trapz(time,(q_temp-mq_c).^2)));
        end
    else
        for k = 1:N
            q_c = q(:,k); mq_c = mu(:,i-1);
            if (option.alignment)
                gamt = optimum_reparam(mq_c,q_c,time,lambda,option.method, ...
                    mu_f(1,i-1), f(1,k));
            else
                gamt = linspace(0,1,M);
            end
            f_temp = warp_f_gamma(f(:,k),gamt,time);
            q_temp = f_to_srvf(f_temp,time);
            Dy(i,k) = sqrt(sum(trapz(time,(q_temp-mq_c).^2)));
        end
    end
    idxold = idx;
    idx = find(rand < cumsum(Dy(i,:)/sum(Dy(i,:))),1);
    Dy(i,idxold)=max(Dy(i,:));
    mu(:,i) = q(:,idx);
    mu_f(:,i) = f(:,idx);
end
lambda = max(min(Dy(2:end,:),[],1));

end

