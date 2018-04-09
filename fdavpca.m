classdef fdavpca
    %fdavpca A class to provide vertical fPCA
    % -------------------------------------------------------------------------
    % This class provides vertical fPCA using the
    % SRVF framework
    %
    % Usage:  obj = fdavpca(warp_data)
    %
    % where:
    %   warp_data - fdawarp object of aligned data
    %
    %
    % fdavpca Properties:
    %   warp_data - fdawarp class with alignment data
    %   q_pca - srvf principal directions
    %   f_pca - f principal directions
    %   latent - latent values
    %   coef - prinicapl coefficients
    %   id - point used for f(0)
    %   mqn - mean srvf
    %   U - eigenvectors
    %   stds - geodesic directions
    %
    %
    % fdavpca Methods:
    %   fdavpca - class constructor
    %   calc_fpca - perform vertical fPCA
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        warp_data % fdawarp class with alignment data
        q_pca     % srvf principal directions
        f_pca     % f principal directions
        latent    % latent values
        coef      % prinicapl coefficients
        id        % point used for f(0)
        mqn       % mean srvf
        U         % eigenvectors
        stds      % geodesic directions
    end
    
    methods
        function obj = fdavpca(fdawarp)
            %fdavpca Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
        end
        
        function obj = calc_fpca(obj,no,id)
            % CALC_FPCA Vertical Functional Principal Component Analysis
            % -------------------------------------------------------------------------
            % This function calculates vertical functional principal component analysis
            % on aligned data
            %
            % Usage: obj.calc_fpca(no)
            %        obj.calc_fpca(no,id)
            %
            % Inputs:
            % warp_data: struct from time_warping of aligned data
            % no: number of principal components to extract
            % id: point to use for f(0) (default = midpoint)
            %
            % Output:
            % fdavpca object
            fn = obj.warp_data.fn;
            t = obj.warp_data.time;
            qn = obj.warp_data.qn;
            
            if nargin < 2
                no = 3;
                id = round(length(t)/2);
                obj.id = id;
            elseif nargin < 3
                id = round(length(t)/2);
                obj.id = id;
            else
                obj.id = id;
            end
            
            % Parameters
            coefs = -2:2;
            NP = 1:no;  % number of principal components
            Nstd = length(coefs);
            
            % FPCA
            mq_new = mean(qn,2);
            m_new = sign(fn(obj.id,:)).*sqrt(abs(fn(obj.id,:)));  % scaled version
            obj.mqn = [mq_new; mean(m_new)];
            K = cov([qn;m_new]');
            
            [obj.U,S,~] = svd(K);
            s = diag(S);
            stdS = sqrt(s);
            s = s(1:no);
            obj.U = obj.U(:,1:no);
            
            % compute the PCA in the q domain
            obj.q_pca = zeros(length(mq_new)+1,Nstd,no);
            for k = NP
                for i = 1:Nstd
                    obj.q_pca(:,i,k) = [mq_new; mean(m_new)] + coefs(i)*stdS(k)*obj.U(:,k);
                end
            end
            
            % compute the correspondence to the original function domain
            obj.f_pca = zeros(length(mq_new),Nstd,no);
            for k = NP
                for i = 1:Nstd
                    obj.f_pca(:,i,k) = cumtrapzmid(t,obj.q_pca(1:end-1,i,k).*abs(obj.q_pca(1:end-1,i,k)),sign(obj.q_pca(end,i,k)).*(obj.q_pca(end,i,k).^2), obj.id);
                end
                fbar = mean(fn,2);
                fsbar = mean(obj.f_pca(:,:,k),2);
                err = repmat(fbar-fsbar,1,Nstd);
                obj.f_pca(:,:,k) = obj.f_pca(:,:,k) + err;
            end
            
            % coefficients
            c = zeros(size(qn,2),no);
            for jj = NP
                for ii = 1:size(qn,2)
                    c(ii,jj) = dot([qn(:,ii);m_new(ii)]-obj.mqn,obj.U(:,jj));
                end
            end
            
            obj.latent = s;
            obj.coef = c;
            obj.stds = coefs;
        end
        
        function plot(obj)
            % plot plot elastic vertical fPCA results
            % -------------------------------------------------------------
            % Usage: obj.plot()
            cl = 'rbgmc';
            [~, ~, p1] = size(obj.q_pca);
            num_plot = ceil(p1/3);
            k = 1;
            for ii = 1:num_plot
                if (k > size(obj.q_pca,3))
                    break
                end
                figure;
                for k1 = 1:3
                    k = k1+(ii-1)*3;
                    subplot(2,3,k1);
                    if (k > size(obj.q_pca,3))
                        break
                    end
                    for i = 1:length(obj.stds)
                        plot(obj.warp_data.time, obj.q_pca(1:end-1,i,k), cl(i), 'linewidth', 2); hold on;
                    end
                    title(['q domain: PD ' num2str(k)], 'fontsize', 14);
                    subplot(2,3,k1+3);
                    for i = 1:length(obj.stds)
                        plot(obj.warp_data.time, obj.f_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
                    end
                    title(['f domain: PD ' num2str(k)], 'fontsize', 14);
                end
            end
            cumm_coef = 100*cumsum(obj.latent)./sum(obj.latent);
            figure
            plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
        end
    end
end
