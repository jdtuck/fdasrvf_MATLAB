classdef fdacurve
    %fdacurve A class to provide registration of curves in R^n using SRVF
    % ---------------------------------------------------------------------
    % This class provides alignment methods for curves in R^n using the
    % SRVF framework
    
    properties
        beta            % (n,T,K) matrix defining n dimensional curve on T samples with K curves
        q               % (n,T,K) matrix defining n dimensional srvf on T samples with K srvfs
        betan           % aligned curves
        qn              % aligned srvfs
        basis           % calculated basis
        beta_mean       % karcher mean curve
        q_mean          % karcher mean srvf
        gams            % warping functions
        v               % shooting vectors
        C               % karcher covariance
        s               % pca singular values
        U               % pca singular vectors
        pca             % principal directions
        coef            % pca coefficients
        closed          % closed curve if true
        qun             % cost function
        samples         % random samples
        gamr            % random warping functions
        cent            % center
        scale           % scale
        E               % energy
        len             % length of curves
        len_q           % length of SRVFs
        mean_scale      % mean length
        mean_scale_q    % mean length SRVF
        center          % centering of curves done
    end
    
    methods
        function obj = fdacurve(beta, closed, N, scale, center)
            %fdacurve Construct an instance of this class
            % Input:
            %   beta: (n,T,K) matrix defining n dimensional curve on T samples with K curves
            %   closed: true or false if closed curve
            %   N: resample curve to N points (default = T)
            %   scale: include scale (true/false (default))
            %   center: center curve (true (default)/false)
            
            if nargin < 3
                N = size(beta,2);
                scale = false;
                center = true;
            elseif nargin < 4
                scale = false;
                center = true;
            elseif nargin < 5
                center = true;
            end
            obj.scale = scale;
            obj.closed = closed;
            
            K = size(beta,3);
            n = size(beta,1);
            q = zeros(n,N,K);
            beta1 = zeros(n,N,K);
            cent1 = zeros(n,K);
            len1 = zeros(1,K);
            lenq1 = zeros(1,K);
            for ii = 1:K
                if size(beta,2) ~= N
                    beta1(:,:,ii) = ReSampleCurve(beta(:,:,ii),N,closed);
                else
                    beta1(:,:,ii) = beta(:,:,ii);
                end
                if (center)
                    a=-calculateCentroid(beta1(:,:,ii));
                else
                    a=zeros(n,1);
                end
                beta1(:,:,ii) = beta1(:,:,ii) + repmat(a,1,N) ;
                [q(:,:,ii), len1(ii), lenq1(ii)] = curve_to_q(beta1(:,:,ii),closed);
                cent1(:,ii) = -a;
            end
            obj.q = q;
            obj.beta = beta1;
            obj.cent = cent1;
            obj.len = len1;
            obj.len_q = lenq1;
            obj.center = center;
        end
        
        function obj = karcher_mean(obj,option)
            % KARCHER_MEAN Calculate karcher mean of group of curves
            % -------------------------------------------------------------
            % This function aligns a collection of functions using the
            % square-root velocity framework
            %
            % Usage:  obj.karcher_mean()
            %         obj.karcher_mean(option)
            %
            % Input:
            %
            % default options
            % option.reparam = true; % computes optimal reparamertization
            % option.rotation = true; % computes optimal rotation
            % option.parallel = 0; % turns on MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.MaxItr = 20;  % maximum iterations
            % option.method = 'DP'; % reparam method
            % controls which optimization method (default="DP") options are
            % Dynamic Programming ("DP") and Riemannian BFGS
            % ("RBFGSM")
            %
            % Output:
            % fdacurve object
            
            if nargin < 2
                option.reparam = true;
                option.rotation = true;
                option.parallel = 0;
                option.closepool = 0;
                option.MaxItr = 20;
                option.method = 'DP';
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
            
            % Initialize mu as one of the shapes
            shape=1;
            mu=obj.q(:,:,shape);
            iter = 1;
            T = size(mu,2);
            n = size(mu,1);
            K = size(obj.q,3);
            gamma = zeros(T,K);
            sumd = zeros(1,option.MaxItr+1);
            normvbar = zeros(1,option.MaxItr+1);
            v1 = zeros(size(obj.q));
            tolv=10^-4;
            told=5*10^-3;
            delta=0.5;
            
            % Compute the Karcher mean
            fprintf('Computing Karcher mean of %d curves in SRVF space...\n',K);
            while iter<=option.MaxItr
                fprintf('updating step: r=%d\n', iter);
                if iter == option.MaxItr
                    fprintf('maximal number of iterations is reached. \n');
                end
                mu=mu/sqrt(InnerProd_Q(mu,mu));
                if obj.closed
                    obj.basis=findBasisNormal(mu);
                end
                
                sumv=zeros(n,T);
                sumd(1) = Inf;
                sumd(iter+1)=0;
                sumnd_t = 0;
                if option.parallel
                    parfor i=1:K
                        q1=obj.q(:,:,i);
                        
                        % Compute shooting vector from mu to q_i
                        [qn_t,~,gamI] = Find_Rotation_and_Seed_unique(mu,q1,option.reparam,option.rotation,obj.closed,option.method);
                        qn_t = qn_t/sqrt(InnerProd_Q(qn_t,qn_t));
                        
                        gamma(:,i) = gamI;
                        
                        q1dotq2=InnerProd_Q(mu,qn_t);
                        
                        % Compute shooting vector
                        if q1dotq2>1
                            q1dotq2=1;
                        end
                        
                        d = acos(q1dotq2);
                        
                        u=qn_t-q1dotq2*mu;
                        normu=sqrt(InnerProd_Q(u,u));
                        if normu>10^-4
                            w=u*acos(q1dotq2)/normu;
                        else
                            w=zeros(size(qn_t));
                        end
                        
                        % Project to tangent space of manifold to obtain v_i
                        if obj.closed
                            v1(:,:,i)=projectTangent(w,q1,obj.basis);
                        else
                            v1(:,:,i)=w;
                        end
                        
                        sumv=sumv+v1(:,:,i);
                        sumnd_t=sumnd_t+d^2;
                    end
                else
                    for i=1:K
                        q1=obj.q(:,:,i);
                        
                        % Compute shooting vector from mu to q_i
                        [qn_t,~,gamI] = Find_Rotation_and_Seed_unique(mu,q1,option.reparam,option.rotation,obj.closed,option.method);
                        qn_t = qn_t/sqrt(InnerProd_Q(qn_t,qn_t));
                        
                        gamma(:,i) = gamI;
                        
                        q1dotq2=InnerProd_Q(mu,qn_t);
                        
                        % Compute shooting vector
                        if q1dotq2>1
                            q1dotq2=1;
                        end
                        
                        d = acos(q1dotq2);
                        
                        u=qn_t-q1dotq2*mu;
                        normu=sqrt(InnerProd_Q(u,u));
                        if normu>10^-4
                            w=u*acos(q1dotq2)/normu;
                        else
                            w=zeros(size(qn_t));
                        end
                        
                        % Project to tangent space of manifold to obtain v_i
                        if obj.closed
                            v1(:,:,i)=projectTangent(w,q1,obj.basis);
                        else
                            v1(:,:,i)=w;
                        end
                        
                        sumv=sumv+v1(:,:,i);
                        sumnd_t=sumnd_t+d^2;
                    end
                end
                sumd(iter+1) = sumnd_t;
                
                % Compute average direction of tangent vectors v_i
                vbar=sumv/K;
                normvbar(iter)=sqrt(InnerProd_Q(vbar,vbar));
                normv=normvbar(iter);
                
                if (sumd(iter)-sumd(iter+1)) < 0
                    break
                elseif (normv>tolv) && abs(sumd(iter+1)-sumd(iter))>told
                    % Update mu in direction of vbar
                    mu=cos(delta*normvbar(iter))*mu+sin(delta*normvbar(iter))*vbar/normvbar(iter);
                    
                    % Project the updated mean to the affine (or closed) shape manifold
                    if obj.closed
                        mu = ProjectC(mu);
                    end
                    
                    x=q_to_curve(mu);
                    if (obj.center)
                        a=-calculateCentroid(x);
                        betamean=x+repmat(a,1,T);
                    else
                        betamean = x;
                    end
                
                else
                    break
                end
                
                iter=iter+1;
                
            end

            % compute average length
            if obj.scale
                % compute geometric mean
                obj.mean_scale = (prod(obj.len))^(1/length(obj.len));
                obj.mean_scale_q = (prod(obj.len_q))^(1/length(obj.len_q));
                betamean = obj.mean_scale.*betamean; 
            end

            % align to mean
            betan1 = obj.beta;
            qn1 = obj.q;
            if option.parallel
                parfor i=1:K
                    q1=obj.q(:,:,i);
                    beta1 = betan1(:,:,i);
                    
                    % Compute shooting vector from mu to q_i
                    [~,R,gamI] = Find_Rotation_and_Seed_unique(mu,q1,option.reparam,option.rotation,obj.closed,option.method);
                    beta1 = R*beta1;
                    beta1n = warp_curve_gamma(beta1,gamI);
                    q1n = curve_to_q(beta1n);
                    
                    % Find optimal rotation
                    [qn1(:,:,i),R] = Find_Best_Rotation(mu,q1n);
                    betan1(:,:,i) = R*beta1n;
                end
                obj.betan = betan1;
                obj.qn = qn1;
            else
                obj.betan = obj.beta;
                obj.qn = obj.q;
                for i=1:K
                    q1=obj.q(:,:,i);
                    beta1 = obj.beta(:,:,i);
                    
                    % Compute shooting vector from mu to q_i
                    [~,R,gamI] = Find_Rotation_and_Seed_unique(mu,q1,option.reparam,option.rotation,obj.closed,option.method);
                    beta1 = R*beta1;
                    beta1n = warp_curve_gamma(beta1,gamI);
                    q1n = curve_to_q(beta1n);
                    
                    % Find optimal rotation
                    [obj.qn(:,:,i),R] = Find_Best_Rotation(mu,q1n);
                    obj.betan(:,:,i) = R*beta1n;
                end
            end
            
            if option.parallel == 1 && option.closepool == 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end
            
            obj.beta_mean = betamean;
            obj.q_mean = mu;
            obj.gams = gamma;
            obj.v = v1;
            obj.qun = sumd(1:iter);
            obj.E = normvbar(1:iter-1);
            
        end
        
        function obj = karcher_cov(obj)
            % KARCHER_COV Calculate karcher covariance
            % -------------------------------------------------------------
            % This function aligns a collection of functions using the
            % square-root velocity framework
            %
            % Usage:  obj.karcher_mean()
            %
            % Input:
            %
            % Output:
            % fdacurve object
            
            [M,N,K] = size(obj.v);
            if obj.scale
                tmpv = zeros(M*N+1,K);
            else
                tmpv = zeros(M*N,K);
            end
            for i = 1:K
                tmp = obj.v(:,:,i);
                if obj.scale
                    tmpv(:,i) = [tmp(:); obj.len(i)];
                else
                    tmpv(:,i) = tmp(:);
                end
            end
            
            if obj.scale
                VM = mean(obj.v,3);
                VM = [VM(:); obj.mean_scale];
                tmpv = tmpv - repmat(VM,1,size(tmpv,2));
                obj.C = cov(tmpv');
            else
                obj.C = cov(tmpv');
            end
            
        end
        
        function obj = shape_pca(obj, no)
            % SHAPE_PCA Calculate Principal Component Analysis
            % -------------------------------------------------------------
            % This function aligns a collection of functions using the
            % square-root velocity framework
            %
            % Usage:  obj.shape_pca()
            %         obj.shape_pca(10)
            %
            % Input:
            % no: number of principal components (default = 10)
            %
            % Output:
            % fdacurve object
            if nargin < 2
                no = 10;
            end
            
            if isempty(obj.C)
                obj = karcher_cov(obj);
            end
            
            % svd
            [U1,S,~] = svd(obj.C);
            ss = diag(S);
            obj.U = U1(:,1:no);
            obj.s = ss(1:no);
            
            % express shapes as coefficients
            K = size(obj.beta,3);
            if obj.scale
                VM = mean(obj.v,3);
                VM = [VM(:); obj.mean_scale];
            else
                VM = mean(obj.v,3);
                VM = VM(:);
            end
            x = zeros(no,K);
            for ii = 1:K
                if obj.scale
                    tmpv = obj.v(:,:,ii);
                    tmpv = [tmpv(:); obj.len(ii)];
                else
                    tmpv = obj.v(:,:,ii);
                    tmpv = tmpv(:);
                end
                x(:,ii) = obj.U'*(tmpv - VM);
            end
            obj.coef = x;
            
            % principal modes of variability
            VM = mean(obj.v,3);
            if obj.scale
                VM = [VM(:); obj.mean_scale];
            else
                VM = VM(:);
            end
            [n,T,~] = size(obj.beta);
            p = zeros(n,T,no,10);
            for j = 1:no
                for i=1:10
                    tmp = VM + 0.5*(i-5)*sqrt(obj.s(j))*obj.U(:,j);
                    [m,n] = size(obj.q_mean);
                    if obj.scale
                        tmp_scale = tmp(end);
                        tmp = tmp(1:end-1);
                    else
                        tmp_scale = 1;
                    end
                    v1 = reshape(tmp,m,n);
                    q2n = ElasticShooting(obj.q_mean,v1);
                    
                    p(:,:,j,i) = q_to_curve(q2n,tmp_scale);
                end
            end
            obj.pca = p;
        end
        
        function obj = sample_shapes(obj, N, m)
            % SAMPLE_SHAPES Sample shapes from generative model
            % -------------------------------------------------------------
            % This function samples shapes from a generative model using a
            % wrapped normal density on shape pca space
            %
            % Usage:  obj.sample_shapes()
            %
            % Input:
            %   N: number of sample shapes
            %   m: number of principal components
            %
            % Output:
            % fdacurve object
            
            if nargin < 3
                N = 10;
                m = 3;
            elseif nargin < 2
                N = 10;
            end
            
            % calculate lengths
            len = zeros(1,size(obj.beta,3));
            T = size(obj.beta,2);
            for ii=1:size(obj.beta,3)
                p = obj.beta(:,:,ii);
                v1=gradient(p,1/(T-1));
                len(ii) = sum(sqrt(sum(v1.*v1)))/T;
            end
            
            if (length(unique(len))==1)
                scale = 1;
            else
                pd = makedist('uniform',min(len),max(len));
            end
            
            % random warping
            gam_s = randomGamma(obj.gams.',N);
            
            [U, S, ~] = svd(obj.C);
            
            % number of shapes in calculating exp mapping
            if obj.closed
                n1 = 10;
            else
                n1 = 2;
            end
            epsilon=1/(n1-1);
            
            n = size(obj.beta_mean,1);
            T = size(obj.beta_mean,2);
            
            obj.samples = zeros(n,T,N);
            for i = 1:N
                v1 = zeros(n,T);
                for dir = 1:m
                    Utemp = reshape(U(:,dir),T,n').';
                    v1 = v1 + randn * sqrt(S(dir,dir))*Utemp;
                end
                
                % q = exp_mu(v) using n1 steps
                q1 = obj.q_mean;
                for j = 1:n1-1
                    normv = sqrt(InnerProd_Q(v1,v1));
                    
                    if normv < 1e-4
                        q2 = mu;
                    else
                        q2 = cos(epsilon*normv)*q1+sin(epsilon*normv)*v1/normv;
                        if obj.closed
                            q2 = ProjectC(q2);
                        end
                    end
                    
                    % parallel translate tangent vector
                    v1=v1-2*InnerProd_Q(v1,q2)/InnerProd_Q(q1+q2,q1+q2)*(q1+q2);
                    
                    q1 = q2;
                end
                
                % random scale
                if (length(unique(len))>1)
                    scale = random(pd);
                end
                beta1s = q_to_curve(q2,scale);
                a = -calculateCentroid(beta1s);
                tmp_beta = beta1s + repmat(a,1,T);
                
                obj.samples(:,:,i) = warp_curve_gamma(tmp_beta,gam_s(i,:))+repmat(mean(obj.cent,2),1,T);
                obj.gamr = gam_s;
                
            end
            
            
            
        end
        
        function plot(obj, color)
            % plot plot curve mean results
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            if nargin < 2
                color = false;
            end
            figure(1);clf;hold all;
            [n,T,K] = size(obj.beta);
            for ii = 1:K
                if n == 2
                    plot(obj.beta(1,:,ii),obj.beta(2,:,ii))
                elseif n == 3
                    if color
                        plot3var(obj.beta(1,:,ii),obj.beta(2,:,ii),obj.beta(3,:,ii),linspace(0,1,T),1)
                    else
                        plot3(obj.beta(1,:,ii),obj.beta(2,:,ii),obj.beta(3,:,ii))
                    end
                else
                    error('Can''t plot dimension > 3')
                end
                title('Curves')
            end
            axis equal ij off;
            
            if (~isempty(obj.gams))
                figure(2); clf
                if n == 2
                    plot(obj.beta_mean(1,:),obj.beta_mean(2,:))
                elseif n == 3
                    plot3(obj.beta_mean(1,:),obj.beta_mean(2,:),obj.beta_mean(3,:))
                else
                    error('Can''t plot dimension > 3')
                end
                title('Karcher Mean')
                axis equal ij off;
                
                figure(3); clf;
                M = size(obj.beta,2);
                plot((0:M-1)/(M-1), obj.gams, 'linewidth', 1);
                axis square;
                grid on;
                title('Warping functions', 'fontsize', 16);
            end
            
            if (~isempty(obj.samples))
                figure(4);clf;hold all;
                K = size(obj.samples,3);
                n = size(obj.samples,1);
                for ii = 1:K
                    if n == 2
                        plot(obj.samples(1,:,ii),obj.samples(2,:,ii))
                    elseif n == 3
                        if color
                            plot3var(obj.samples(1,:,ii),obj.samples(2,:,ii),obj.samples(3,:,ii),linspace(0,1,T),1)
                        else
                            plot3(obj.samples(1,:,ii),obj.samples(2,:,ii),obj.samples(3,:,ii))
                        end
                    else
                        error('Can''t plot dimension > 3')
                    end
                    title('Sample Curves')
                end
                axis equal ij off;
            end
        end
        
        function plot_pca(obj, n)
            
            if nargin < 2
                n = 4;
            end
            if isempty(obj.s)
                error('Calculate PCA')
            end
            figure(5)
            plot(cumsum(obj.s)/sum(obj.s)*100)
            title('Variability Explained')
            xlabel('PC')
            
            % plot principal modes of variability
            VM = mean(obj.v,3);
            if obj.scale
                VM = [VM(:); obj.mean_scale];
            else
                VM = VM(:);
            end
            for j = 1:n
                figure(20+j); clf; hold on;
                for i=1:10
                    tmp = VM + 0.5*(i-5)*sqrt(obj.s(j))*obj.U(:,j);
                    [m,n] = size(obj.q_mean);
                    if obj.scale
                        tmp_scale = tmp(end);
                        tmp = tmp(1:end-1);
                    else
                        tmp_scale = 1;
                    end
                    v1 = reshape(tmp,m,n);
                    q2n = ElasticShooting(obj.q_mean,v1);
                    
                    p = q_to_curve(q2n,tmp_scale);
                    if obj.scale
                        mv = 0.2*obj.mean_scale;
                    else
                        mv = 0.2;
                    end
                    if i == 5
                        plot(mv*i + p(1,:),p(2,:), 'k','LineWidth',3);
                    else
                        plot(mv*i + p(1,:),p(2,:),'LineWidth',2);
                    end
                    axis equal ij off;
                end
                title(sprintf('PC: %d',j))
            end
        end
        
    end
end

function plot3var(x,y,z,c,lwd)
surface('XData',[x(:) x(:)], 'YData',[y(:) y(:)], 'ZData',[z(:) z(:)], ...
    'CData',[c(:) c(:)], 'FaceColor','none', 'EdgeColor','interp', ...
    'Marker','none','linew',lwd)
end

function plot_curve2(beta)
plot(beta(1,:),beta(2,:),'LineWidth',1.5);
axis equal ij off;
end