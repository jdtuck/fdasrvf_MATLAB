classdef fdavpca
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        warp_data
        q_pca
        f_pca
        latent
        coef
        id
        mqn
        U
        stds
    end
    
    methods
        function obj = fdavpca(fdawarp)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.warp_data = fdawarp;
        end
        
        function obj = calc_fpca(obj,no,id)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fn = obj.warp_data.fn;
            t = obj.warp_data.time;
            qn = obj.warp_data.qn;
            
            if nargin < 3
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

