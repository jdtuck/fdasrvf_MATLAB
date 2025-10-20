classdef fdahpns
    %fdahpns A class to provide horizontal fPNS
    % -------------------------------------------------------------------------
    % This class provides horizontal fPNS using the
    % SRVF framework
    %
    % Usage:  obj = fdahpns(warp_data)
    %
    % where:
    %   warp_data - fdawarp object of aligned data
    %
    %
    % fdahpns Properties:
    %   warp_data - fdawarp class with alignment data
    %   gam_pns - warping functions principal directions
    %   psi_pns - srvf principal directions
    %   latent - latent values
    %   coef - coeficients
    %   stds - principal directions
    %   new_coef - principal coefficients of new data
    %
    %
    % fdahpns Methods:
    %   fdahpns - class constructor
    %   calc_fpns - perform horizontal fPCA
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  13-Sep-2025

    properties
        warp_data % fdawarp class with alignment data
        gam_pns   % warping functions principal directions
        psi_pns   % srvf principal directions
        latent    % latent values
        pc        % eigenvectors
        coef      % coeficients
        cumvar    % cummulative variance
        no        % no of principal spheres
        PNS       % PNS structure
        stds       % principal directions
        new_coef  % principal coefficients of new data 
    end

    methods
        function obj = fdahpns(fdawarp)
            %fdahpca Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
        end

        function obj = calc_fpns(obj,no,stds)
            % calc_fpns Horizontal Functional Principal Component Analysis
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
            % stds: number of standard deviations along geodesic to compute
            %       (default = -2:2)
            %
            % Outputs:
            % fdahpca object
            arguments
                obj
                no = 3;
                stds = -2:2;
            end
            gam = obj.warp_data.gam;
            n = size(gam, 2);
            M = size(gam, 1);

            [resmat, obj.PNS] = PNS_warping(gam);

            % proportion of variance explained
            varPNS = sum(abs(resmat.^2), 2) / n;
            cumvarPNS = cumsum(varPNS);
            propcumPNS = cumvarPNS / cumvarPNS(end);

            % Parameters
            obj.stds = stds;

            % TFPCA
            obj.psi_pns = zeros(length(obj.stds),M,no);
            obj.gam_pns = zeros(length(obj.stds),M,no);
            for j=1:no      % three components
                std1 = std(resmat(j, :));
                mean1 = mean(resmat(j, :));
                dirtmp = obj.stds * std1 + mean1;
                restmp = zeros(size(resmat,1), length(obj.stds));
                restmp(j, :) = dirtmp;
                PCvec = fastPNSe2s(restmp, stds);
                obj.psi_pns(:, :, j) = (PCvec * stds.radius);
                for k=1:length(obj.stds)
                    gam0 = cumtrapz(linspace(0,1,size(gam,1)),obj.psi_pns(k,:,j).*obj.psi_pns(k,:,j));
                    obj.gam_pns(k,:,j) = (gam0-gam0(1))/(gam0(end)-gam0(1));
                end
            end

            obj.cumvar = propcumPNS;
            obj.no = no;
            obj.coef = resmat;
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

            psi = zeros(M, n);
            time = linspace(0,1,M);
            binsize = mean(diff(time));
            for i = 1:n
                psi(:, i) = sqrt(gradietn(gam(:, i), binsize));
            end

            pnsdat = psi./repmat(sqrt(sum(psi.^2)),d,1);

            dat = obj.PNS.basisu' * pnsdat;
            
            resmat = fastPNSe2s(dat, obj.PNS);

            obj.new_coef = resmat;
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
            [~, T, p1] = size(obj.gam_pns);
            num_plot = ceil(p1/3);
            j = 1;
            for ii = 1:num_plot
                if (j > size(obj.gam_pns,3))
                    break
                end
                figure;
                for j1 = 1:3
                    j = j1 + (ii-1)*3;
                    if (j > size(obj.gam_pns,3))
                        break
                    end
                    subplot(1,3,j1);
                    for k = 1:length(obj.stds)
                        plot(linspace(0,1,T), obj.gam_pns(k,:,j), cl(k), 'linewidth', 2); hold on;
                    end
                    plot(linspace(0,1,T), obj.gam_pns(idx,:,j), 'k', 'linewidth', 2)
                    axis([0 1 0 1]);
                    title(['PD ' num2str(j)], 'fontsize', 14);
                end
            end

            figure
            plot(obj.cumvar, 'Color', cl(1));title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
        end
    end
end
