classdef ampbox
    %ampbox Produce amplitude box plots
    % -------------------------------------------------------------------------
    % This class provides amplitude boxplot for functional data using the
    % SRVF framework
    %
    % Usage:  obj = ampbox(warp_data)
    %
    % where:
    %   warp_data - fdawarp class of aligned data
    %
    %
    % ampbox Properties:
    %   warp_data - fdawarp class with alignment data
    %   Q1 - First quartile
    %   Q3 - Second quartile
    %   Q1a - First quantile based on alpha
    %   Q3a - Second quantile based on alpha
    %   minn - minimum extreme function
    %   maxx - maximum extreme function
    %   outlier_index - indexes of outlier functions
    %   f_median - median function
    %   q_median - median srvf
    %   plt - surface plot mesh
    %
    %
    % ampbox Methods:
    %   ampbox - class constructor
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
        f_median       % median function
        q_median       % median srvf
        plt            % surface plot mesh
        Q1_index       % index of quartiles
        Q3_index       % index of quartiles
        Q1a_index      % index of quantiles
        Q3a_index      % index of quantiles
        dist           % distances
    end
    
    methods
        function obj = ampbox(fdawarp)
            %ampbox Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping_median!');
            end
            
            if (~strcmpi(fdawarp.type,'median'))
                error('Please align fdawarp class using time_warping_median!');
            end
            obj.warp_data = fdawarp;
        end
        
        function obj = construct_boxplot(obj, alpha, k_a)
            % CONSTRUCT_BOXPLOT constructs the amplitude boxplot
            % -------------------------------------------------------------------------
            %
            % Usage:  obj.construct_boxplot(alpha, k_a)
            %
            % Input:
            % alpha: quantile value (e.g.,=.05, i.e., 95\%)
            % ka: scalar for outlier cutoff (e.g.,=1)
            %
            % Output: structure containing
            % ampbox object
            
            if obj.warp_data.rsamps
                f_tilde = obj.warp_data.fs;
                obj.f_median = obj.warp_data.fmean;
                q_tilde = obj.warp_data.qs;
                obj.q_median = obj.warp_data.mqn;
                t = obj.warp_data.time;
            else
                f_tilde = obj.warp_data.fn;
                obj.f_median = obj.warp_data.fmean;
                q_tilde = obj.warp_data.qn;
                obj.q_median = obj.warp_data.mqn;
                t = obj.warp_data.time;
            end
            
            [~, N] = size(f_tilde);
            lambda = 0.5;
                        
            % compute amplitude distances
            dy = zeros(1,N);
            for i = 1:N
                dy(i) = sqrt(trapz(t,(obj.q_median-q_tilde(:,i)).^2));
            end
            obj.dist = dy;
            [~, dy_ordering] = sort(dy);
            CR_50 = dy_ordering(1:ceil(N/2));       % 50% Central Region
            m = max(dy(CR_50));                     % Maximal amplitude distance within 50% Central Region
            
            % identify amplitude quartiles
            angle = zeros(length(CR_50), length(CR_50));
            energy = zeros(length(CR_50), length(CR_50));
            for i = 1:(length(CR_50)-1)
                for j = (i+1):length(CR_50)
                    q1 = q_tilde(:,CR_50(i)) - obj.q_median;
                    q3 = q_tilde(:,CR_50(j)) - obj.q_median;
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda) * (dy(CR_50(i))/m + dy(CR_50(j))/m) - lambda * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            obj.Q1_index = CR_50(maxloc_row);
            obj.Q3_index = CR_50(maxloc_col);
            Q1_q = q_tilde(:,obj.Q1_index);
            Q3_q = q_tilde(:,obj.Q3_index);
            obj.Q1 = f_tilde(:,obj.Q1_index);
            obj.Q3 = f_tilde(:,obj.Q3_index);
            
            % identify amplitude quantiles
            [~, dy_ordering] = sort(dy);
            CR_alpha = dy_ordering(1:round(N*(1-alpha)));       % (1-alpha)% Central Region
            m = max(dy(CR_alpha));                     % Maximal amplitude distance within (1-alpha)% Central Region
            angle = zeros(length(CR_alpha), length(CR_alpha));
            energy = zeros(length(CR_alpha), length(CR_alpha));
            for i = 1:(length(CR_alpha)-1)
                for j = (i+1):length(CR_alpha)
                    q1 = q_tilde(:,CR_alpha(i)) - obj.q_median;
                    q3 = q_tilde(:,CR_alpha(j)) - obj.q_median;
                    q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
                    q3=q3/sqrt(trapz(t,q3.^2));
                    angle(i,j)=trapz(t,q1.*q3);
                    energy(i,j) = (1-lambda) * (dy(CR_alpha(i))/m + dy(CR_alpha(j))/m) - lambda * (angle(i,j) + 1);
                end
            end
            [~, maxloc] = max(energy(:));
            [maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);
            
            obj.Q1a_index = CR_alpha(maxloc_row);
            obj.Q3a_index = CR_alpha(maxloc_col);
            Q1a_q = q_tilde(:,obj.Q1a_index);
            Q3a_q = q_tilde(:,obj.Q3a_index);
            obj.Q1a = f_tilde(:,obj.Q1a_index);
            obj.Q3a = f_tilde(:,obj.Q3a_index);
            
            % compute amplitude whiskers
            IQR = dy(obj.Q1_index)+dy(obj.Q3_index);
            v1 = Q1_q - obj.q_median;
            v3 = Q3_q - obj.q_median;
            upper_q = Q3_q + k_a * IQR * v3 / sqrt(trapz(t,v3.^2));
            lower_q = Q1_q + k_a * IQR * v1 / sqrt(trapz(t,v1.^2));
            
            upper_dis = sqrt(trapz(t,(upper_q - obj.q_median).^2));
            lower_dis = sqrt(trapz(t,(lower_q - obj.q_median).^2));
            whisker_dis = max([lower_dis upper_dis]);
            
            % identify amplitude outliers
            obj.outlier_index = [];
            for i = 1:N
                if dy(dy_ordering(N+1-i)) > whisker_dis
                    obj.outlier_index = [obj.outlier_index; dy_ordering(N+1-i)];
                else
                    break
                end
            end
            
            % identify amplitude extremes
            distance_to_upper=inf(1,N);
            distance_to_lower=inf(1,N);
            out_50_CR = setdiff(setdiff((1:N), CR_50), obj.outlier_index);
            for i = 1:length(out_50_CR)
                j = out_50_CR(i);
                distance_to_upper(j) = sqrt(trapz(t,(upper_q - q_tilde(:,j)).^2));
                distance_to_lower(j) = sqrt(trapz(t,(lower_q - q_tilde(:,j)).^2));
            end
            [~, max_index] = min(distance_to_upper);
            [~, min_index] = min(distance_to_lower);
            min_q = q_tilde(:,min_index);
            max_q = q_tilde(:,max_index);
            obj.minn = f_tilde(:,min_index);
            obj.maxx = f_tilde(:,max_index);
            
            s = linspace(0,1,100);
            Fs2 = zeros(length(t), 595);
            Fs2(:,1) = (1-s(1)) * obj.minn + s(1) * obj.Q1;    % Final surface plot
            for j=2:100
                Fs2(:,j) = (1-s(j)) * obj.minn + s(j) * obj.Q1a;
                Fs2(:,99+j) = (1-s(j)) * obj.Q1a + s(j) * obj.Q1;
                Fs2(:,198+j) = (1-s(j)) * obj.Q1 + s(j) * obj.f_median;
                Fs2(:,297+j) = (1-s(j)) * obj.f_median + s(j) * obj.Q3;
                Fs2(:,396+j) = (1-s(j)) * obj.Q3 + s(j) * obj.Q3a;
                Fs2(:,495+j) = (1-s(j)) * obj.Q3a + s(j) * obj.maxx;
            end
            d1=sqrt(trapz(t,(obj.q_median-Q1_q).^2));
            d1a=sqrt(trapz(t,(Q1_q-Q1a_q).^2));
            dl=sqrt(trapz(t,(Q1a_q-min_q).^2));
            d3=sqrt(trapz(t,(obj.q_median-Q3_q).^2));
            d3a=sqrt(trapz(t,(Q3_q-Q3a_q).^2));
            du=sqrt(trapz(t,(Q3a_q-max_q).^2));
            part1=linspace(-d1-d1a-dl,-d1-d1a,100);
            part2=linspace(-d1-d1a,-d1,100);
            part3=linspace(-d1,0,100);
            part4=linspace(0,d3,100);
            part5=linspace(d3,d3+d3a,100);
            part6=linspace(d3+d3a,d3+d3a+du,100);
            allparts=[part1,part2(2:100),part3(2:100),part4(2:100),part5(2:100),part6(2:100)];
            [U,V]=meshgrid(t,allparts);
            U=U';
            V=V';
            
            obj.plt.U=U;
            obj.plt.V=V;
            obj.plt.Fs2 = Fs2;
            obj.plt.allparts = allparts;
            obj.plt.d1 = d1;
            obj.plt.d1a = d1a;
            obj.plt.dl = dl;
            obj.plt.d3 = d3;
            obj.plt.d3a = d3a;
            obj.plt.du = du;
            obj.plt.Q1q = Q1a_q;
            obj.plt.Q3q = Q3a_q;
            
        end
        
        function plot(obj)
            % plot plot box plot and surface plot
            % -------------------------------------------------------------------------
            % Usage: obj.plot()
            figure(310); clf;
            t = obj.warp_data.time;
            M = length(t);
            plot(t, obj.warp_data.fmean, 'black','linewidth', 2);
            hold on;
            plot(t, obj.Q1, 'blue','linewidth', 2);
            plot(t, obj.Q3, 'blue', 'linewidth', 2);
            plot(t, obj.Q1a, 'green', 'linewidth', 2);
            plot(t, obj.Q3a, 'green', 'linewidth', 2);
            plot(t, obj.minn,'red', 'linewidth',2);
            plot(t, obj.maxx,'red', 'linewidth',2);
            xlim([t(1) t(end)]);
            ylim auto;
            
            
            figure(311); clf;
            surf(obj.plt.U,obj.plt.V,obj.plt.Fs2);
            hold on;
            shading flat;
            plot3(t,zeros(1,M),obj.f_median,'k','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1,M,1),obj.Q1,'b','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1-obj.plt.d1a,M,1),obj.Q1a,'g','LineWidth',3)
            plot3(t,repmat(-obj.plt.d1-obj.plt.d1a-obj.plt.dl,M,1),obj.minn,'r','LineWidth',3)
            plot3(t,repmat(obj.plt.d3,M,1),obj.Q3,'b','LineWidth',3)
            plot3(t,repmat(obj.plt.d3+obj.plt.d3a,M,1),obj.Q3a,'g','LineWidth',3)
            plot3(t,repmat(obj.plt.d3+obj.plt.d3a+obj.plt.du,M,1),obj.maxx,'r','LineWidth',3)
        end
    end
end

