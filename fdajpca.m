classdef fdajpca
    %fdajpca A class to provide a joint fPCA
    % -------------------------------------------------------------------------
    % This class provides joint fPCA using the
    % SRVF framework
    %
    % Usage:  obj = fdajpca(warp_data)
    %
    % where:
    %   warp_data - fdawarp object of aligned data
    %
    %
    % fdajpca Properties:
    %   warp_data - fdawarp class with alignment data
    %   q_pca - srvf principal directions
    %   f_pca - f principal directions
    %   latent - latent values
    %   coef - prinicapl coefficients
    %   id - point used for f(0)
    %   mqn - mean srvf
    %   U - eigenvectors
    %   mu_psi - mean psi
    %   mu_g - mean g
    %   C - scaling value
    %   stds - geodesic directions
    %
    %
    % fdajpca Methods:
    %   fdajpca - class constructor
    %   calc_fpca - perform vertical fPCA
    %   plot - plot results and functions in object
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  18-Mar-2018
    
    properties
        warp_data % fdawarp class with alignment data
        q_pca     % srvf principal directions
        f_pca     % f principal directions
        latent    % latent values
        coef      % prinicapl coefficients
        id        % point used for f(0)
        mqn       % mean srvf
        U         % eigenvectors
        mu_psi    % mean psi
        mu_g      % mean g
        C         % scaling value
        stds      % geodesic directions
    end
    
    methods
        function obj = fdajpca(fdawarp)
            %fdajpca Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
        end
        
        function obj = calc_fpca(obj,no,id)
            % CALC_FPCA Joint Functional Principal Component Analysis
            % -------------------------------------------------------------------------
            % This function calculates joint functional principal component analysis
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
            % fdajpca object
            fn = obj.warp_data.fn;
            time = obj.warp_data.time;
            qn = obj.warp_data.qn;
            q0 = obj.warp_data.q0;
            gam = obj.warp_data.gam;
            if nargin < 2
                no = 3;
                id = round(length(time)/2);
                obj.id = id;
            elseif nargin < 3
                id = round(length(time)/2);
                obj.id = id;
            else
                obj.id = id;
            end
            
            [M, ~] = size(qn);
            
            % set up for fPCA in q-space
            mq_new = mean(qn,2);
            id = round(length(time)/2);
            m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));  % scaled version
            obj.mqn = [mq_new; mean(m_new)];
            qn1 = [qn; m_new];
            
            % calculate vector space of warping functions
            [obj.mu_psi,~,~,vec] = SqrtMean(gam.');
            vec = vec.';
            
            % joint fPCA
            f1 = @(x)find_C(x, qn1, vec, q0, no, obj.mu_psi);
            obj.C = fminbnd(f1, 0, 1e4);
            
            [~, ~, a, obj.U, s, obj.mu_g] = jointfPCAd(qn1, vec, obj.C, no, obj.mu_psi);
            
            % geodesic paths
            ci = [-1,0,1];
            obj.stds = ci;
            obj.q_pca = zeros(M, length(ci), no);
            obj.f_pca = zeros(M, length(ci), no);
            
            for j = 1:no
                for i = 1:length(ci)
                    qhat = obj.mqn + obj.U(1:(M+1),j) * ci(i)*sqrt(s(j));
                    vechat = obj.U((M+2):end,j) * (ci(i)*sqrt(s(j)))/obj.C;
                    psihat = exp_map(obj.mu_psi,vechat);
                    gamhat = cumtrapz(linspace(0,1,M), psihat.*psihat);
                    gamhat = (gamhat - min(gamhat))/(max(gamhat)-min(gamhat));
                    if (sum(vechat)==0)
                        gamhat = linspace(0,1,M);
                    end
                    
                    fhat = cumtrapzmid(time, qhat(1:M).*abs(qhat(1:M)),sign(qhat(M+1)).*(qhat(M+1)^2), id);
                    obj.f_pca(:,i,j) = warp_f_gamma(fhat, gamhat, linspace(0,1,M));
                    obj.q_pca(:,i,j) = warp_q_gamma(qhat(1:M), gamhat, linspace(0,1,M));
                end
            end
            
            obj.coef = a;
            obj.latent = s;
        end
        
        function plot(obj)
            % plot plot elastic vertical fPCA results
            % -------------------------------------------------------------
            % Usage: obj.plot()
            cl = 'rbgmc';
            [~, ~, p1] = size(obj.q_pca);
            time = obj.warp_data.time;
            ci = obj.stds;
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
                    for i = 1:length(ci)
                        plot(time, obj.q_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
                    end
                    title(['q domain: PD ' num2str(k)], 'fontsize', 14);
                    subplot(2,3,k1+3);
                    for i = 1:length(ci)
                        plot(time, obj.f_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
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

function [qhat, gamhat, a, U, s, mu_g] = jointfPCAd(qn, vec, C, m, mu_psi)
[M, N] = size(qn);
g = [qn; C*vec];
mu_q = mean(qn,2);
mu_g = mean(g,2);

K = cov(g.');

[U,S,~] = svd(K);
s = diag(S);

a = zeros(N,m);
for i = 1:N
    for j = 1:m
        a(i,j) = (g(:,i)-mu_g).'*U(:,j);
    end
end

qhat = repmat(mu_q,1,N) + U(1:M,1:m) * a.';
vechat = U((M+1):end,1:m) * (a.'/C);
psihat = zeros(M-1,N);
gamhat = zeros(M-1,N);
for ii = 1:N
    psihat(:,ii) = exp_map(mu_psi,vechat(:,ii));
    gam_tmp = cumtrapz(linspace(0,1,M-1), psihat(:,ii).*psihat(:,ii));
    gamhat(:,ii) = (gam_tmp - min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
end

U = U(:,1:m);
s = s(1:m);

end

function [out] = find_C(C, qn, vec, q0, m, mu_psi)
[qhat, gamhat, ~, ~, ~, ~] = jointfPCAd(qn, vec, C, m, mu_psi);
[M, N] = size(qn);
time = linspace(0,1,M-1);

d = zeros(1,N);
for i = 1:N
    tmp = warp_q_gamma(qhat(1:(M-1),i), invertGamma(gamhat(:,i)), time);
    d(i) = sum(trapz(time,(tmp-q0(:,i)).^2));
end

out = sum(d.^2)/N;

end

