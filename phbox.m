classdef phbox
    %phbox Construct phase box plots
    % -------------------------------------------------------------------------
    % This class provides phase boxplot for functional data using the
    % SRVF framework
    %
    % Usage:  obj = phbox(warp_data)
    %
    % where:
    %   warp_data - fdawarp class of aligned data
    %
    %
    % phbox Properties:
    %   warp_data      % fdawarp class with alignment data
    %   Q1             % First quartile
    %   Q3             % Second quartile
    %   Q1a            % First quantile based on alpha
    %   Q3a            % Second quantile based on alpha
    %   minn           % minimum extreme function
    %   maxx           % maximum extreme function
    %   outlier_index  % indexes of outlier functions
    %   median_x       % median warping function
    %   psi_median     % median srvf of warping function
    %   plt            % surface plot mesh
    %
    %
    % phbox Methods:
    %   phbox - class constructor
    %   construct_boxplot - construct boxplot
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  15-Mar-2018
    
    properties
        warp_data      % fdawarp class with alignment data
        Q1             % First quartile
        Q3             % Second quartile
        Q1a            % First quantile based on alpha
        Q3a            % Second quantile based on alpha
        minn           % minimum extreme function
        maxx           % maximum extreme function
        outlier_index  % indexes of outlier functions
        median_x       % median warping function
        psi_median     % median srvf of warping function
        plt            % surface plot mesh
        Q1_index       % index of quartiles
        Q3_index       % index of quartiles
        Q1a_index      % index of quantiles
        Q3a_index      % index of quantiles
        dist           % distances
    end
    
    methods
        function obj = phbox(fdawarp)
            %phbox Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            if (~strcmpi(fdawarp.type,'median'))
                error('Please align fdawarp class using time_warping_median!');
            end
            obj.warp_data = fdawarp;
        end
        
        function obj = construct_boxplot(obj, alpha, k_p)
            % CONSTRUCT_BOXPLOT constructs the phase boxplot
            % -------------------------------------------------------------------------
            %
            % Usage:  obj.construct_boxplot(k_p, alpha)
            %
            % Input:
            % alpha quantile: value (e.g.,=.05, i.e., 95\%)
            % k_p: scalar for outlier cutoff (e.g.,=1)
            %
            % Output: structure containing
            % phbox object
            
            if obj.warp_data.rsamps
                gam = obj.warp_data.gams;
            else
                gam = obj.warp_data.gam;
            end
            
            [M, N] = size(gam);
            t = linspace(0,1,M);
            lambda = 0.5;
            
            % compute phase median
            [obj.median_x, obj.psi_median, psi] = SqrtMedian(gam);
            
            % compute phase distances
            binsize = mean(diff(t));
            dx = zeros(1,N);
            v = zeros(M,N);
            for i = 1:N
                psi(:,i) = sqrt(gradient(gam(:,i),binsize));
                v(:,i) = inv_exp_map(obj.psi_median,psi(:,i));
                dx(i) = sqrt(trapz(t,v(:,i).^2));
            end
            obj.dist = dx;
            [~, dx_ordering] = sort(dx);
            CR_50 = dx_ordering(1:ceil(N/2));   % 50% Central Region
            m = max(dx(CR_50));                 % Maximal phase distance within 50% Central Region
            
            % identify phase quartiles
            angle = zeros(length(CR_50), length(CR_50));
            energy = zeros(length(CR_50), length(CR_50));
            for i = 1:(length(CR_50)-1)
                for j = (i+1):length(CR_50)
                    q1 = v(:,CR_50(i));
                    q3 = v(:,CR_50(j));
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda) * (dx(CR_50(i))/m + dx(CR_50(j))/m) - lambda * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            obj.Q1_index = CR_50(maxloc_row);
            obj.Q3_index = CR_50(maxloc_col);
            obj.Q1 = gam(:, obj.Q1_index);
            obj.Q3 = gam(:, obj.Q3_index);
            Q1_psi = sqrt(gradient(obj.Q1,1/(M-1)))';
            Q3_psi = sqrt(gradient(obj.Q3,1/(M-1)))';
            
            % identify phase quantiles
            [~, dx_ordering] = sort(dx);
            CR_alpha = dx_ordering(1:round(N*(1-alpha)));   % 50% Central Region
            m = max(dx(CR_alpha));                 % Maximal phase distance within 50% Central Regionles
            
            angle = zeros(length(CR_alpha), length(CR_alpha));
            energy = zeros(length(CR_alpha), length(CR_alpha));
            for i = 1:(length(CR_alpha)-1)
                for j = (i+1):length(CR_alpha)
                    q1 = v(:,CR_alpha(i));
                    q3 = v(:,CR_alpha(j));
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda) * (dx(CR_alpha(i))/m + dx(CR_alpha(j))/m) - lambda * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            obj.Q1a_index = CR_alpha(maxloc_row);
            obj.Q3a_index = CR_alpha(maxloc_col);
            obj.Q1a = gam(:, obj.Q1a_index);
            obj.Q3a = gam(:, obj.Q3a_index);
            Q1a_psi = sqrt(gradient(obj.Q1a,1/(M-1)))';
            Q3a_psi = sqrt(gradient(obj.Q3a,1/(M-1)))';
            
            % check quartile and quatnile going same direction
            tst = trapz(t, v(:,obj.Q1a_index).*v(:,obj.Q1_index));
            if (tst < 0)
                obj.Q1a = gam(:,obj.Q3a_index);
                obj.Q3a = gam(:,obj.Q1a_index);
            end
            
            % compute phase whiskers
            IQR = dx(obj.Q1_index) + dx(obj.Q3_index);
            v3 = v(:, obj.Q3_index);
            v1 = v(:, obj.Q1_index);
            upper_v = v3 + k_p *IQR * v3 / sqrt(trapz(t,v3.^2));
            lower_v = v1 + k_p *IQR * v1 / sqrt(trapz(t,v1.^2));
            
            whisker_dis = max(lower_v,upper_v);
            
            % identify phase outliers
            obj.outlier_index = [];
            for i = 1:N
                if dx(dx_ordering(N+1-i)) > whisker_dis
                    obj.outlier_index = [obj.outlier_index; dx_ordering(N+1-i)];
                else
                    break
                end
            end
            
            % identify phase extremes
            distance_to_upper=inf(1,N);
            distance_to_lower=inf(1,N);
            out_50_CR = setdiff(setdiff(1:N, CR_50), obj.outlier_index);
            for i = 1:length(out_50_CR)
                j = out_50_CR(i);
                distance_to_upper(j) = sqrt(trapz(t,(upper_v - v(:,j)).^2));
                distance_to_lower(j) = sqrt(trapz(t,(lower_v - v(:,j)).^2));
            end
            [~, max_index] = min(distance_to_upper);
            [~, min_index] = min(distance_to_lower);
            obj.minn = gam(:,min_index);
            obj.maxx = gam(:,max_index);
            min_psi = psi(:,min_index);
            max_psi = psi(:,max_index);
            
            s = linspace(0,1,100);
            t = t(:);
            obj.median_x = obj.median_x(:);
            Fs2 = zeros(length(t), 595);
            Fs2(:,1) = (1-s(1)) * (obj.minn-t) + s(1) * (obj.Q1-t);
            for j=2:100
                Fs2(:,j) = (1-s(j)) * (obj.minn-t) + s(j) * (obj.Q1a-t);
                Fs2(:,99+j) = (1-s(j)) * (obj.Q1a-t) + s(j) * (obj.Q1-t);
                Fs2(:,198+j) = (1-s(j)) * (obj.Q1-t) + s(j) * (obj.median_x-t);
                Fs2(:,297+j) = (1-s(j)) * (obj.median_x-t) + s(j) * (obj.Q3-t);
                Fs2(:,396+j) = (1-s(j)) * (obj.Q3-t) + s(j) * (obj.Q3a-t);
                Fs2(:,495+j) = (1-s(j)) * (obj.Q3a-t) + s(j) * (obj.maxx-t);
            end
            Q1_psi = Q1_psi(:);
            Q1a_psi = Q1a_psi(:);
            Q3_psi = Q3_psi(:);
            Q3a_psi = Q3a_psi(:);
            d1=acos(trapz(t,obj.psi_median.*Q1_psi));
            d1a=acos(trapz(t,Q1_psi.*Q1a_psi));
            dl=acos(trapz(t,Q1a_psi.*min_psi));
            d3=acos(trapz(t,obj.psi_median.*Q3_psi));
            d3a=acos(trapz(t,Q3_psi.*Q3a_psi));
            du=acos(trapz(t,Q3a_psi.*max_psi));
            part1=linspace(-d1-d1a-dl,-d1-d1a,100);
            part2=linspace(-d1-d1a,-d1,100);
            part3=linspace(-d1,0,100);
            part4=linspace(0,d3,100);
            part5=linspace(d3,d3+d3a,100);
            part6=linspace(d3+d3a,d3+d3a+du,100);
            allparts=[part1,part2(2:100),part3(2:100),part4(2:100),part5(2:100),part6(2:100)];
            [U,V]=meshgrid(linspace(0,1,M),allparts);
            U=U.';
            V=V.';
            
            obj.plt.U=U;
            obj.plt.V=V;
            obj.plt.allparts = allparts;
            obj.plt.Fs2 = Fs2;
            obj.plt.d1 = d1;
            obj.plt.d1a = d1a;
            obj.plt.dl = dl;
            obj.plt.d3 = d3;
            obj.plt.d3a = d3a;
            obj.plt.du = du;
            obj.plt.Q1_psi = Q1a_psi;
            obj.plt.Q3_psi = Q3a_psi;
            
        end
        
        function plot(obj)
            % plot plot box plot and surface plot
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            [M, ~] = size(obj.warp_data.gam);
            t = linspace(0,1,M).';
            figure(410); clf;
            plot(t, obj.median_x, 'black','linewidth', 2);
            hold on;
            plot(t, obj.Q1, 'blue','linewidth', 2);
            plot(t, obj.Q3, 'blue', 'linewidth', 2);
            plot(t, obj.Q1a, 'green', 'linewidth', 2);
            plot(t, obj.Q3a, 'green', 'linewidth', 2);
            plot(t, obj.maxx,'red','linewidth',2);
            plot(t, obj.minn,'red','linewidth',2);
            axis square;
            axis([0,1,0,1]);
            
            figure(416); clf;
            surf(obj.plt.U,obj.plt.V,obj.plt.Fs2);
            hold on;
            shading flat;
            plot3(t,zeros(1,M),obj.median_x - t,'k','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1,M,1),obj.Q1 - t,'b','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1-obj.plt.d1a,M,1),obj.Q1a - t,'g','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1-obj.plt.d1a-obj.plt.dl,M,1),obj.minn -t,'r','LineWidth',3)
            plot3(t,repmat(obj.plt.d3,M,1),obj.Q3 - t,'b','LineWidth',3)
            plot3(t,repmat(obj.plt.d3+obj.plt.d3a,M,1),obj.Q3a - t,'g','LineWidth',3)
            plot3(t,repmat(obj.plt.d3+obj.plt.d3a+obj.plt.du,M,1),obj.maxx - t,'r','LineWidth',3)
            axis square;
        end
    end
end

