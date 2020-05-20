classdef fdacurve
    %fdacurve A class to provide registration of curves in R^n using SRVF
    % -------------------------------------------------------------------------
    % This class provides alignment methods for curves in R^n using the
    % SRVF framework

    properties
        beta      % (n,T,K) matrix defining n dimensional curve on T samples with K curves
        q         % (n,T,K) matrix defining n dimensional srvf on T samples with K srvfs
        betan     % aligned curves
        qn        % aligned srvfs
        basis     % calculated basis
        beta_mean % karcher mean curve
        q_mean    % karcher mean srvf
        gams      % warping functions
        v         % shooting vectors
        C         % karcher covariance
        closed    % closed curve if true
        qun       % cost function
        samples   % random samples
        gamr      % random warping functions
        cent      % center
        scale     % scale
    end

    methods
        function obj = fdacurve(beta, closed, N, scale)
            %fdacurve Construct an instance of this class
            % Input:
            %   beta: (n,T,K) matrix defining n dimensional curve on T samples with K curves
            %   closed: true or false if closed curve
            %   N: resample curve to N points
            %   scale: scale curve to length 1 (true/false)

            if nargin < 4
                scale = true;
            end
            obj.scale = scale;
            obj.closed = closed;

            K = size(beta,3);
            n = size(beta,1);
            q = zeros(n,N,K);
            beta1 = zeros(n,N,K);
            cent1 = zeros(n,K);
            for ii = 1:K
                beta1(:,:,ii) = ReSampleCurve(beta(:,:,ii),N,closed);
                a=-calculateCentroid(beta1(:,:,ii));
                beta1(:,:,ii) = beta1(:,:,ii) + repmat(a,1,N) ;
                q(:,:,ii) = curve_to_q(beta1(:,:,ii),obj.scale,closed);
                cent1(:,ii) = -a;
            end
            obj.q = q;
            obj.beta = beta1;
            obj.cent = cent1;
        end

        function obj = karcher_mean(obj,option)
            % KARCHER_MEAN Calculate karcher mean of group of curves
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the
            % square-root velocity framework
            %
            % Usage:  obj.karcher_mean()
            %         obj.karcher_mean(option)
            %
            % Input:
            %
            % default options
            % option.parallel = 0; % turns on MATLAB parallel processing (need
            % parallel processing toolbox)
            % option.closepool = 0; % determines wether to close matlabpool
            % option.MaxItr = 20;  % maximum iterations
            %
            % Output:
            % fdacurve object

            if nargin < 2
                option.parallel = 0;
                option.closepool = 0;
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
                        [qn_t,~,gamI] = Find_Rotation_and_Seed_unique(mu,q1,true,obj.closed);
                        if obj.scale
                            qn_t = qn_t/sqrt(InnerProd_Q(qn_t,qn_t));
                        end
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
                        [qn_t,~,gamI] = Find_Rotation_and_Seed_unique(mu,q1,true,obj.closed);
                        if obj.scale
                            qn_t = qn_t/sqrt(InnerProd_Q(qn_t,qn_t));
                        end
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

                if (normv>tolv) && abs(sumd(iter+1)-sumd(iter))>told
                    % Update mu in direction of vbar
                    mu=cos(delta*normvbar(iter))*mu+sin(delta*normvbar(iter))*vbar/normvbar(iter);

                    % Project the updated mean to the affine (or closed) shape manifold
                    if obj.closed
                        mu = ProjectC(mu);
                    end

                    x=q_to_curve(mu);
                    a=-calculateCentroid(x);
                    betamean=x+repmat(a,1,T);

                else
                    break
                end

                iter=iter+1;

            end

            betan1 = obj.beta;
            qn1 = obj.q;
            if option.parallel
                parfor i=1:K
                    q1=obj.q(:,:,i);
                    beta1 = betan1(:,:,i);

                    % Compute shooting vector from mu to q_i
                    [~,R,gamI] = Find_Rotation_and_Seed_unique(mu,q1,true,obj.closed);
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
                    [~,R,gamI] = Find_Rotation_and_Seed_unique(mu,q1,true,obj.closed);
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
            obj.E = normvbar(iter);

        end

        function obj = karcher_cov(obj)
            % KARCHER_COV Calculate karcher covariance
            % -------------------------------------------------------------------------
            % This function aligns a collection of functions using the
            % square-root velocity framework
            %
            % Usage:  obj.karcher_mean()
            %
            % Input:
            %
            % Output:
            % fdacurve object

            [n,T,K] = size(obj.v);

            obj.C=zeros(n*T);
            for i=1:K
                w=obj.v(:,:,i);
                [m1,n1] = size(w);
                w = reshape(w',1,m1*n1);
                obj.C=obj.C+w'*w;
            end
            obj.C=obj.C/(K-1);
        end

        function obj = sample_shapes(obj, N, m)
            % SAMPLE_SHAPES Sample shapes from generative model
            % -------------------------------------------------------------------------
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
            end
        end
    end
end

function plot3var(x,y,z,c,lwd)
surface('XData',[x(:) x(:)], 'YData',[y(:) y(:)], 'ZData',[z(:) z(:)], ...
    'CData',[c(:) c(:)], 'FaceColor','none', 'EdgeColor','interp', ...
    'Marker','none','linew',lwd)
end
