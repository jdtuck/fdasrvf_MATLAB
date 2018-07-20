classdef fdacurve
    %fdacurve A class to provide registration of curves in R^n using SRVF
    % -------------------------------------------------------------------------
    % This class provides alignment methods for curves in R^n using the
    % SRVF framework
    
    properties
        beta      % (n,T,K) matrix defining n dimensional curve on T samples with K curves
        q         % (n,T,K) matrix defining n dimensional srvf on T samples with K srvfs
        basis     % calculated basis
        beta_mean % karcher mean curve
        q_mean    % karcher mean srvf
        gams      % warping functions
        v         % shooting vectors
        C         % karcher covariance
        closed    % closed curve if true
    end
    
    methods
        function obj = fdacurve(beta, closed, N)
            %fdacurve Construct an instance of this class
            % Input:
            %   beta: (n,T,K) matrix defining n dimensional curve on T samples with K curves
            %   closed: true or false if closed curve
            %   N: resample curve to N points
            
            obj.closed = closed;
            
            K = size(beta,3);
            n = size(beta,1);
            q = zeros(n,N,K);
            beta1 = zeros(n,N,K);
            for ii = 1:K
                beta1(:,:,ii) = ReSampleCurve(beta(:,:,ii),N);
                q(:,:,ii) = curve_to_q(beta1(:,:,ii));
            end
            obj.q = q;
            obj.beta = beta1;
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
            K = size(obj.q,3);
            gamma = zeros(T,K);
            sumd = zeros(1,option.MaxItr+1);
            normvbar = zeros(1,option.MaxItr+1);
            v1 = zeros(size(obj.q));
            tolv=10^-4;
            told=5*10^-3;
            delta=0.5;
            
            % Compute the Karcher mean
            while iter<=option.MaxItr
                mu=mu/sqrt(InnerProd_Q(mu,mu));
                if obj.closed
                    obj.basis=Basis_Normal_A(mu);
                end
                
                sumv=zeros(2,T);
                sumd(iter+1)=0;
                
                for i=1:K
                    q1=obj.q(:,:,i);
                    
                    % Compute shooting vector from mu to q_i
                    [qn,~,gamI] = Find_Rotation_and_Seed_unique(mu,q1,true,obj.closed);
                    gamma(:,i) = gamI;
                    [qn,~] = Find_Best_Rotation(mu,qn);
                    
                    q1dotq2=InnerProd_Q(mu,qn);
                    
                    % Compute shooting vector
                    if q1dotq2>1
                        q1dotq2=1;
                    end
                    
                    d = acos(q1dotq2);
                    
                    u=qn-q1dotq2*mu;
                    normu=sqrt(InnerProd_Q(u,u));
                    if normu>10^-4
                        w=u*acos(q1dotq2)/normu;
                    else
                        w=zeros(2,T);
                    end
                    
                    % Project to tangent space of manifold to obtain v_i
                    if obj.closed
                        v1(:,:,i)=projectTangent(w,q1);
                    else
                        v1(:,:,i)=w;
                        
                    end
                    
                    sumv=sumv+v1(:,:,i);
                    sumd(iter+1)=sumd(iter+1)+d^2;
                end
                
                % Compute average direction of tangent vectors v_i
                vbar=sumv/T;
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
            
            obj.beta_mean = betamean;
            obj.q_mean = mu;
            obj.gams = gamma;
            obj.v = v1;
            
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
            
            [~,T,K] = size(obj.v);
            
            obj.C=zeros(2*T);
            for i=1:K
                w=obj.v(:,:,i);
                w=[w(1,:) w(2,:)];
                obj.C=obj.C+w'*w;
            end
            obj.C=obj.C/(K-1);
        end
        
        function plot(obj)
            % plot plot curve mean results
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            figure(1);clf;hold all;
            K = size(obj.beta,3);
            for ii = 1:K
                plot(obj.beta(1,:,ii),obj.beta(2,:,ii))
                title('Curves')
            end
            
            if (~isempty(obj.gams))
                figure(2); clf
                plot(obj.beta_mean(1,:),obj.beta_mean(2,:))
                title('Karcher Mean')
                
                figure(3); clf;
                M = size(obj.beta,2);
                plot((0:M-1)/(M-1), obj.gams, 'linewidth', 1);
                axis square;
                title('Warping functions', 'fontsize', 16);
            end
        end
    end
end

