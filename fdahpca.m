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
    %   vm - mean of shooting vectors
    %   stds - principal directions
    %   new_coef - principal coefficients of new data  
    %   log_der - log-derivative transform used
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
        vm        % mean of shooting vectors
        stds      % principal directions
        new_coef  % principal coefficients of new data 
        log_der   % log-derivative transform used
    end
    
    methods
        function obj = fdahpca(fdawarp, log_der)
            %fdahpca Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            arguments
                fdawarp
                log_der = false
            end
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
            obj.log_der = log_der;
        end
        
        function obj = calc_fpca(obj, no, var_exp, stds)
            % calc_fpca Horizontal Functional Principal Component Analysis
            % -------------------------------------------------------------------------
            % This function calculates vertical functional principal component analysis
            % on aligned data
            %
            % Usage: obj.calc_fpca(no)
            %
            % Inputs:
            % no: number of principal components to extract
            % var_exp: compute no based on value percent variance explained
            %          (example: 0.95)
            % stds: number of standard deviations along geodesic to compute
            %       (default = -2:2)
            %
            % Outputs:
            % fdahpca object
            
            arguments
                obj
                no = 3;
                var_exp = NaN;
                stds = -2:2;
            end
            gam = obj.warp_data.gam;

            idx = find(stds == 0);
            if isempty(idx)
                error("stds needs to contain 0")
            end

            M = length(obj.warp_data.time);
            if ~isnan(var_exp)
                if var_exp > 1
                    error("var_exp is greater than 1")
                end
                no = M;
            end

            if obj.log_der
                obj.vec = gam_to_h(gam);
                obj.mu = mean(obj.vec, 2);
            else
                [obj.mu,~,obj.vec] = SqrtMean(gam);
            end

            K = cov(obj.vec');
            
            [obj.U,S,~] = svd(K);
            Sig = diag(S);
            if ~isnan(var_exp)
                cumm_coef = cumsum(Sig) / sum(Sig);
                no = find(cumm_coef >= var_exp, 'first');
            end
            T = size(obj.vec,1);
            obj.U = obj.U(:,1:no);
            Sig = Sig(1:no);
            
            % Parameters
            obj.stds = stds;
            
            % TFPCA
            obj.psi_pca = zeros(length(obj.stds),length(obj.mu),no);
            obj.gam_pca = zeros(length(obj.stds),length(obj.mu),no);
            v = zeros(5,length(obj.mu),3);
            for j=1:no      % three components
                for k=1:length(obj.stds)   % -2, -1, 0, 1, 2 std from the mean
                    v(k,:,j) = obj.stds(k)*sqrt(Sig(j))*obj.U(:,j)';
                    if obj.log_der
                        obj.psi_pca(k,:,j) = v(k,:,j);
                        gam0 = h_to_gam(obj.psi_pca(k,:,j));
                    else
                        vn = norm(v(k,:,j))/sqrt(T);
                        if vn < 0.0001
                            obj.psi_pca(k,:,j) = obj.mu;
                        else
                            obj.psi_pca(k,:,j) = cos(vn)*obj.mu.' + sin(vn)*v(k,:,j)/vn;
                        end
                        gam0 = cumtrapz(linspace(0,1,size(gam,1)),obj.psi_pca(k,:,j).*obj.psi_pca(k,:,j));
                    end
                    
                    obj.gam_pca(k,:,j) = (gam0-gam0(1))/(gam0(end)-gam0(1));
                end
            end
            
            % coefficients
            c = zeros(size(gam,1),no);
            obj.vm = mean(obj.vec,2);
            for jj = 1:no
                for ii = 1:size(gam,2)
                    c(ii,jj) = (obj.vec(:,ii)-obj.vm).'*obj.U(:,jj);
                end
            end
            
            obj.latent = Sig;
            obj.coef = c;
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
            gam = zeros(M, n);
            for ii = 1:n
                gam(:, ii) = optimum_reparam(mq, obj.warp_data.time, q1(:, ii));
            end

            no = size(obj.U,2);

            mu_psi = obj.mu;
            vec1 = zeros(M, n);
            psi = zeros(M, n);
            time = linspace(0,1,M);
            binsize = mean(diff(time));
            for i = 1:n
                if obj.log_der
                    out = gam_to_h(gam(:,i));
                else
                    psi(:, i) = sqrt(gradietn(gam(:, i), binsize));
                    [out, ~] = inv_exp_map(mu_psi, psi(:,i));
                end
                vec1(:, i) = out;
            end

            c = zeros(size(gam,1),no);
            for jj = 1:no
                for ii = 1:size(gam,2)
                    c(ii,jj) = (vec1(:,ii)-obj.vm).'*obj.U(:,jj);
                end
            end

            obj.new_coef = c;
        end
        
        function plot(obj)
            % plot plot elastic horizontal fPCA results
            % -------------------------------------------------------------------------
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
            [~, T, p1] = size(obj.gam_pca);
            idx = find(obj.stds == 0);
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
                    for k = 1:length(obj.stds)
                        plot(linspace(0,1,T), obj.gam_pca(k,:,j), 'Color', cl(k), 'linewidth', 2); hold on;
                    end
                    plot(linspace(0,1,T), obj.gam_pca(idx,:,j), 'k', 'linewidth', 2)

                    axis([0 1 0 1]);
                    title(['PD ' num2str(j)], 'fontsize', 14);
                end
            end
            cumm_coef = 100*cumsum(obj.latent)./sum(obj.latent);
            figure
            plot(cumm_coef, 'Color', cl(1));title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
        end
    end
end
