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
    %   new_coef - principal coefficients of new data  
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
        new_coef  % principal coefficients of new data 
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
        
        function obj = calc_fpca(obj,no,id,stds)
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
            % stds: number of standard deviations along geodesic to compute
            %       (default = -2:2)
            %
            % Output:
            % fdavpca object

            arguments
                obj
                no = 3;
                id = round(length(obj.warp_data.time)/2);
                stds = -2:2;
            end
            obj.id = id;
            fn = obj.warp_data.fn;
            t = obj.warp_data.time;
            qn = obj.warp_data.qn;
            
            % Parameters
            NP = 1:no;  % number of principal components
            Nstd = length(stds);
            
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
                    obj.q_pca(:,i,k) = [mq_new; mean(m_new)] + stds(i)*stdS(k)*obj.U(:,k);
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
            obj.stds = stds;
        end

        function obj = project(obj, f)
            % PROJECT Project new data onto fPCA basis
            % -------------------------------------------------------------------------
            % This function project new data onto fPCA basis
            %
            % Usage: obj.project(f)
            %        obj.calc_fpca(no,id)
            %
            % Inputs:
            % f:  array (MxN) of N functions on M time points
            %

            q1 = f_to_srvf(f, obj.warp_data.time);
            M = length(obj.warp_data.time);
            n = size(q1,2);
            mq = obj.warp_data.mqn;
            fn = zeros(M, n);
            qn = zeros(M, n);
            gam = zeros(M, n);
            for ii = 1:n
                gam(:, ii) = optimum_reparam(mq, obj.warp_data.time, q1(:, ii));
                fn(:, ii) = warp_f_gamma(f(:, ii), gam(:, ii), obj.warp_data.time);
                qn(:, ii) = f_to_srvf(fn(:, ii), obj.warp_data.time);
            end

            no = size(obj.U,2);

            m_new = sign(fn(obj.id, :)) .* sqrt(abs(fn(obj.id, :)));

            c = zeros(n,no);
            for jj = 1:no
                for ii = 1:n
                    c(ii,jj) = dot([qn(:,ii);m_new(ii)]-obj.mqn,obj.U(:,jj));
                end
            end

            obj.new_coef = c;
        end
        
        function plot(obj)
            % plot plot elastic vertical fPCA results
            % -------------------------------------------------------------
            % Usage: obj.plot()
            cl = [
                "#66C2A5";
                "#FC8D62";
                "#8DA0CB";
                "#E78AC3";
                "#A6D854";
                "#FFD92F";
                "#E5C494";
                "#B3B3B3"
              ];
            idx = find(obj.stds == 0);
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
                        if i == idx
                            color = [0, 0, 0];
                        else
                            color = cl(i);
                        end
                        plot(obj.warp_data.time, obj.q_pca(1:end-1,i,k), 'Color', color, 'linewidth', 2); hold on;
                    end
                    title(['q domain: PD ' num2str(k)], 'fontsize', 14);
                    subplot(2,3,k1+3);
                    for i = 1:length(obj.stds)
                        if i == idx
                            color = [0, 0, 0];
                        else
                            color = cl(i);
                        end
                        plot(obj.warp_data.time, obj.f_pca(:,i,k), 'Color', color, 'linewidth', 2); hold on;
                    end
                    title(['f domain: PD ' num2str(k)], 'fontsize', 14);
                end
            end
            cumm_coef = 100*cumsum(obj.latent)./sum(obj.latent);
            figure
            plot(cumm_coef, 'Color', cl(1));title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
        end
    end
end
