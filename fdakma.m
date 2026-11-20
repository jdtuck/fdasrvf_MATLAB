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
    %   gamI - inverse Karcher mean warp of each cluster (cell)
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
    %   DPmeans - perform DP-means (number of clusters learned)
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
        gamI         % inverse Karcher mean warp of each cluster (cell)
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
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.alignment = true; % performa alignment
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, DP2, SIMUL, RBFGS)
            % option.MaxItr = 20;  % maximum iterations
            % option.thresh = 0.01; % cost function threshold
            %
            % Output:
            % fdakma object
            
            arguments
                obj
                K
                seeds = [];
                lambda1 = 0;
                option.parallel = 0;
                option.nworkers = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP';
                option.MaxItr = 20;
                option.thresh = 0.01;
            end

            obj.lambda = lambda1;
            
            % time warping on a set of functions
            if option.parallel == 1
                start_parpool(option.nworkers);
            end
            %% Parameters
            
            fprintf('\n lambda = %5.1f \n', obj.lambda);
            
            
            [M, N] = size(obj.f);
            
            if (isempty(seeds))
                template_ind = randsample(N,K);
            else
                template_ind = seeds;
            end
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end

            obj.templates = zeros(M,K);
            for i=1:K
                obj.templates(:,i) = obj.f(:,template_ind(i));
            end
            cluster_id = zeros(1,N);
            obj.qun = zeros(1, option.MaxItr);
            obj.gamI = cell(1,K);

            
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
                    if ~any(id)
                        % empty cluster: nothing to normalize
                        obj.gamI{i} = linspace(0,1,M);
                        continue
                    end
                    ftmp = fn1{i}(:,id);
                    gamtmp = gam1{i}(id,:);
                    % SqrtMeanInverse expects warps as columns (samples x functions)
                    gamIc = SqrtMeanInverse(gamtmp.');
                    obj.gamI{i} = gamIc;
                    N1 = size(ftmp,2);
                    f_temp1 = zeros(M,N1);
                    q_temp1 = zeros(M,N1);
                    gam2 = zeros(N1,M);
                    if option.parallel == 1
                        parfor k = 1:N1
                            f_temp1(:,k) = warp_f_gamma(ftmp(:,k),gamIc,obj.time);
                            q_temp1(:,k) = f_to_srvf(f_temp1(:,k),obj.time);
                            gam2(k,:) = warp_f_gamma(gamtmp(k,:),gamIc,obj.time);
                        end
                    else
                        for k = 1:N1
                            f_temp1(:,k) = warp_f_gamma(ftmp(:,k),gamIc,obj.time);
                            q_temp1(:,k) = f_to_srvf(f_temp1(:,k),obj.time);
                            gam2(k,:) = warp_f_gamma(gamtmp(k,:),gamIc,obj.time);
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
                    if ~any(id)
                        % empty cluster: keep the previous template
                        continue
                    end
                    old_template_q = obj.templates_q(:,i);
                    obj.templates_q(:,i) = mean(qn1{i}(:,id),2);
                    % keep the f-template consistent with the SRVF template
                    obj.templates(:,i) = srvf_to_f(obj.templates_q(:,i),obj.time, ...
                        mean(fn1{i}(1,id)));
                    qun_t(i) = norm(obj.templates_q(:,i)-old_template_q)/norm(old_template_q);
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
                if ~isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            obj.fn = ftmp;
            obj.qn = qtmp;
            obj.gam = gamtmp;
            obj.labels = cluster_id;
            obj.method = option.method;
            obj.qun = obj.qun(1:r);
        end
        
        function obj = DPmeans(obj,lambda1,option)
            % DPmeans Dirichlet Process K-Means clustering and alignment
            % -------------------------------------------------------------------------
            % This function clusters functions and aligns using the elastic square-root
            % slope (srsf) framework using the DP-means algorithm (Kulis and Jordan,
            % 2012), the small-variance limit of a Dirichlet process mixture. A new
            % cluster is opened whenever a function's elastic distance to every
            % existing template exceeds lambda1, so the number of clusters is
            % learned from the data. Each pass assigns and aligns the functions,
            % normalizes the warps in each cluster, and updates the templates as
            % in fdakma.kmeans.
            %
            % Usage:  obj.DPmeans()
            %         obj.DPmeans(lambda1)
            %         obj.DPmeans(lambda1,option)
            %
            % Input:
            % lambda1 : cluster penalty, as an elastic distance (default NaN: chosen
            %           by a k-means++ style heuristic). The objective minimized is
            %           sum of squared distances + lambda1^2 * K
            %
            % default options
            % option.lambda_align = 0; % alignment penalty passed to optimum_reparam
            % option.parallel = 0; % turns offs MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.nworkers = 0; % number of workers for parallel pool (0: prompt)
            % option.alignment = true; % performa alignment
            % option.closepool = 0; % determines wether to close matlabpool
            % option.smooth = 0; % smooth data using standard box filter
            % option.sparam = 25; % number of times to run filter
            % option.method = 'DP'; % optimization method (DP, SIMUL, RBFGS)
            % option.MaxItr = 20;  % maximum iterations
            % option.thresh = 0.01; % relative change in objective for convergence
            %
            % Output:
            % fdakma object

            arguments
                obj
                lambda1 = NaN
                option.lambda_align = 0;
                option.parallel = 0;
                option.nworkers = 0;
                option.alignment = true;
                option.closepool = 0;
                option.smooth = 0;
                option.sparam = 25;
                option.method = 'DP';
                option.MaxItr = 20;
                option.thresh = 0.01;
            end

            % time warping on a set of functions
            if option.parallel == 1
                start_parpool(option.nworkers);
            end

            %% Parameters
            if option.smooth == 1
                obj.f = smooth_data(obj.f, option.sparam);
            end
            [M, N] = size(obj.f);

            %% Compute the q-function of the plot
            q = f_to_srvf(obj.f,obj.time);
            obj.q0 = q;

            % start with one cluster whose template is the function closest
            % to the cross-sectional mean srvf
            mnq = mean(q,2);
            dqq = sqrt(trapz(obj.time,(q - mnq).^2,1));
            [~, min_ind] = min(dqq);
            K = 1;
            obj.templates = obj.f(:,min_ind);
            obj.templates_q = q(:,min_ind);

            obj.lambda = lambda1;
            if isnan(obj.lambda)
                fprintf('\n Initializing Lambda \n');
                Kinit = 4;
                obj.lambda = kpp_init(obj.f,q,obj.time,Kinit,option);
            end

            fprintf('\n lambda = %5.3f \n', obj.lambda);

            cluster_id = ones(1,N);
            dmin = zeros(1,N);
            fa = zeros(M,N);
            qa = zeros(M,N);
            ga = zeros(N,M);
            obj.gamI = {};
            obj.qun = zeros(1,option.MaxItr);

            fprintf('\nInitializing...\n');
            for r = 1:option.MaxItr
                fprintf('updating step: r=%d\n', r);
                if r == option.MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end

                %% Assignment and alignment
                % align every function to every current template
                gam1 = cell(1,K);
                fn1 = cell(1,K);
                qn1 = cell(1,K);
                Dy = zeros(K,N);
                for i = 1:K
                    [gam1{i},fn1{i},qn1{i},Dy(i,:)] = align_to_template(obj.f,q, ...
                        obj.time,obj.templates(:,i),obj.templates_q(:,i),option);
                end

                % sequential DP-means pass: a function farther than lambda from
                % every template opens a new cluster, which later functions in
                % the pass can join
                for k = 1:N
                    [dk, ik] = min(Dy(:,k));
                    if dk > obj.lambda
                        K = K + 1;
                        obj.templates(:,K) = obj.f(:,k);
                        obj.templates_q(:,K) = q(:,k);
                        [gam1{K},fn1{K},qn1{K},Dy(K,:)] = align_to_template(obj.f,q, ...
                            obj.time,obj.templates(:,K),obj.templates_q(:,K),option);
                        ik = K;
                        dk = Dy(K,k);
                    end
                    cluster_id(k) = ik;
                    dmin(k) = dk;
                end

                % drop clusters that lost all their functions
                [used,~,cluster_id] = unique(cluster_id);
                cluster_id = cluster_id(:).';
                K = numel(used);
                obj.templates = obj.templates(:,used);
                obj.templates_q = obj.templates_q(:,used);
                gam1 = gam1(used);
                fn1 = fn1(used);
                qn1 = qn1(used);

                % DP-means objective at this assignment
                obj.qun(r) = sum(dmin.^2) + obj.lambda^2*K;

                % keep each function aligned to its own template
                for k = 1:N
                    ga(k,:) = gam1{cluster_id(k)}(k,:);
                    fa(:,k) = fn1{cluster_id(k)}(:,k);
                    qa(:,k) = qn1{cluster_id(k)}(:,k);
                end

                %% Normalization
                obj.gamI = cell(1,K);
                for i = 1:K
                    id = find(cluster_id == i);
                    gamIc = SqrtMeanInverse(ga(id,:).');
                    obj.gamI{i} = gamIc;
                    for k = id
                        fa(:,k) = warp_f_gamma(fa(:,k),gamIc,obj.time);
                        qa(:,k) = f_to_srvf(fa(:,k),obj.time);
                        ga(k,:) = warp_f_gamma(ga(k,:),gamIc,obj.time);
                    end
                end

                %% Template Identification
                for i = 1:K
                    id = cluster_id == i;
                    obj.templates_q(:,i) = mean(qa(:,id),2);
                    obj.templates(:,i) = srvf_to_f(obj.templates_q(:,i),obj.time, ...
                        mean(fa(1,id)));
                end

                if r > 1 && abs(obj.qun(r)-obj.qun(r-1)) < option.thresh*obj.qun(r-1)
                    break;
                end
            end

            %% Output
            ftmp = cell(1,K);
            qtmp = cell(1,K);
            gamtmp = cell(1,K);
            for i = 1:K
                id = cluster_id == i;
                ftmp{i} = fa(:,id);
                qtmp{i} = qa(:,id);
                gamtmp{i} = ga(id,:);
            end

            if option.parallel == 1 && option.closepool == 1
                if ~isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end

            obj.fn = ftmp;
            obj.qn = qtmp;
            obj.gam = gamtmp;
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
            colors = turbo(K);
            figure(2); clf; hold on
            for k=1:K
                plot(obj.time, obj.templates(:,k), 'Color', colors(k,:));
            end
            title('Cluster Mean Functions');
            
            figure(3); clf;
            plot(obj.time, obj.fn{1}, 'Color', colors(1,:));
            hold on
            for k = 2:K
                plot(obj.time, obj.fn{k}, 'Color', colors(k,:));
            end
            for k = 1:K
                plot(obj.time, obj.templates(:,k), 'Color', colors(k,:));
            end
            title('Clustered Functions');
            
            figure(4); clf;
            plot(obj.time, obj.qn{1}, 'Color', colors(1,:));
            hold on
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

function lambda = kpp_init(f,q,time,K,option)
% KPP_INIT choose the DP-means penalty from a k-means++ seeding
% Picks K centers by k-means++ (sampling proportional to the squared elastic
% distance to the nearest chosen center) and returns the largest distance
% from any function to its nearest center (Kulis and Jordan, 2012).
N = size(q,2);
D = zeros(K,N);
idx = randi(N);
for i = 1:K
    [~,~,~,D(i,:)] = align_to_template(f,q,time,f(:,idx),q(:,idx),option);
    if i < K
        w = min(D(1:i,:),[],1).^2;
        if sum(w) == 0
            D = D(1:i,:);
            break
        end
        idx = find(rand*sum(w) < cumsum(w),1);
    end
end
lambda = max(min(D,[],1));
end

function [gamt,ft,qt,d] = align_to_template(f,q,time,tf,tq,option)
% ALIGN_TO_TEMPLATE align every function to one template
% gamt is N x M, ft and qt are M x N, d is the elastic distance (1 x N)
[M,N] = size(f);
gamt = zeros(N,M);
ft = zeros(M,N);
qt = zeros(M,N);
d = zeros(1,N);
if option.parallel == 1
    nw = Inf;
else
    nw = 0;
end
lam = option.lambda_align;
method = option.method;
alignment = option.alignment;
parfor (k = 1:N, nw)
    if alignment
        g = optimum_reparam(tq,q(:,k),time,lam,method,tf(1),f(1,k));
    else
        g = linspace(0,1,M);
    end
    fw = warp_f_gamma(f(:,k),g,time);
    qw = f_to_srvf(fw,time);
    gamt(k,:) = g;
    ft(:,k) = fw;
    qt(:,k) = qw;
    d(k) = sqrt(trapz(time,(qw-tq).^2));
end
end
