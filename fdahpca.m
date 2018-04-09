classdef fdahpca
    %fdahpca A class to provide horizontal fPCA
    % -------------------------------------------------------------------------
    % This class provides horizontal fPCA using the
    % SRVF framework
    %
    % Usage:  obj = fdahpca(warp_data)
    %
    % where:
    %   warp_data - fdawarp object of aligned data
    %
    %
    % fdahpca Properties:
    %   warp_data - fdawarp class with alignment data
    %   gam_pca - warping functions principal directions
    %   psi_pca - srvf principal directions
    %   latent - latent values
    %   U - eigenvectors
    %   coef - coeficients
    %   vec - shooting vectors
    %   mu - Karcher Mean
    %   tau - principal directions
    %
    %
    % fdahpca Methods:
    %   fdahpca - class constructor
    %   calc_fpca - perform horizontal fPCA
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        warp_data % fdawarp class with alignment data
        gam_pca   % warping functions principal directions
        psi_pca   % srvf principal directions
        latent    % latent values
        U         % eigenvectors
        coef      % coeficients
        vec       % shooting vectors
        mu        % Karcher Mean
        tau       % principal directions
    end
    
    methods
        function obj = fdahpca(fdawarp)
            %fdahpca Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
        end
        
        function obj = calc_fpca(obj,no)
            % calc_fpca Horizontal Functional Principal Component Analysis
            % -------------------------------------------------------------------------
            % This function calculates vertical functional principal component analysis
            % on aligned data
            %
            % Usage: obj.horizFPCA(no)
            %        obj.horizFPCA(no,showplot)
            %
            % Inputs:
            % no: number of principal components to extract
            % showplot: show plots of principal directions (e.g, true)
            %
            % Outputs:
            % fdahpca object
            gam = obj.warp_data.gam;
            if (~exist('no'))
                no = 3;
            end
            
            [obj.mu,~,obj.vec] = SqrtMean(gam.');
            
            K = cov(obj.vec);
            
            [obj.U,S,~] = svd(K);
            Sig = diag(S);
            T = size(obj.vec,2);
            obj.U = obj.U(:,1:no);
            Sig = Sig(1:no);
            
            % Parameters
            obj.tau = 1:5;
            
            % TFPCA
            obj.psi_pca = zeros(length(obj.tau),length(obj.mu),no);
            obj.gam_pca = zeros(length(obj.tau),length(obj.mu),no);
            v = zeros(5,length(obj.mu),3);
            for j=1:no      % three components
                for k=obj.tau   % -2, -1, 0, 1, 2 std from the mean
                    v(k,:,j) = (k-3)*sqrt(Sig(j))*obj.U(:,j)';
                    vn = norm(v(k,:,j))/sqrt(T);
                    if vn < 0.0001
                        obj.psi_pca(k,:,j) = obj.mu;
                    else
                        obj.psi_pca(k,:,j) = cos(vn)*obj.mu.' + sin(vn)*v(k,:,j)/vn;
                    end
                    gam0 = cumtrapz(linspace(0,1,size(gam,1)),obj.psi_pca(k,:,j).*obj.psi_pca(k,:,j));
                    obj.gam_pca(k,:,j) = (gam0-gam0(1))/(gam0(end)-gam0(1));
                end
            end
            
            % coefficients
            c = zeros(size(gam,2),no);
            obj.vec = obj.vec.';
            vm = mean(obj.vec,2);
            for jj = 1:no
                for ii = 1:size(gam,2)
                    c(ii,jj) = (obj.vec(:,ii)-vm).'*obj.U(:,jj);
                end
            end
            
            obj.latent = Sig;
            obj.coef = c;
        end
        
        function plot(obj)
            % plot plot elastic horizontal fPCA results
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            cl = 'rbgmc';
            [~, T, p1] = size(obj.gam_pca);
            num_plot = ceil(p1/3);
            j = 1;
            for ii = 1:num_plot
                if (j > size(obj.gam_pca,3))
                    break
                end
                figure;
                for j1 = 1:3
                    j = j1 + (ii-1)*3;
                    if (j > size(obj.gam_pca,3))
                        break
                    end
                    subplot(1,3,j1);
                    for k = 1:length(obj.tau)
                        plot(linspace(0,1,T), obj.gam_pca(k,:,j), cl(k), 'linewidth', 2); hold on;
                    end
                    axis([0 1 0 1]);
                    title(['PD ' num2str(j)], 'fontsize', 14);
                end
            end
            cumm_coef = 100*cumsum(obj.latent)./sum(obj.latent);
            figure
            plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
        end
    end
end
