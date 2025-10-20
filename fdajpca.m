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
    %   U_q - eigenvectors of function 
    %   U_h - eigenvectors of warping
    %   vec - vector space rep of warping
    %   no  - number of principal coefficients
    %   qn1 - representation of functions in srvf or not 
    %   U1  - eigenvectors of sub
    %   Uz  - eigenvectors of sub
    %   mv  - mean of vector rep
    %   C - scaling value
    %   stds - geodesic directions
    %   new_coef - principal coefficients of new data
    %   new_g - g of new data
    %   srvf - srvf space used
    %   log_der - log-derivative transform used
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
        U_q       % eigenvectors of function 
        U_h       % eigenvectors of warping
        vec       % vector space rep of warping
        no        % number of principal coefficients
        qn1       % representation of functions in srvf or not 
        U1        % eigenvectors of sub
        Uz        % eigenvectors of sub
        mv        % mean of vector rep
        C         % scaling value
        stds      % geodesic directions
        new_coef  % principal coefficients of new data 
        new_h     % vec rep of gamma of new data
        new_qn1   % qn rep of new data
        srvf      % srvf space used
        log_der   % log-derivative transform used
    end
    
    methods
        function obj = fdajpca(fdawarp, log_der)
            %fdajpca Construct an instance of this class
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
        
        function obj = calc_fpca(obj,var_exp,id,stds,srvf)
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
            % var_exp: compute no based on value percent variance explained (default = 0.99)
            % id: point to use for f(0) (default = midpoint)
            % stds: number of standard deviations along geodesic to compute
            %       (default = -2:2)
            % srvf: use srvf (default = T)
            %
            % Output:
            % fdajpca object
            arguments
                obj
                var_exp=0.99;
                id = round(length(obj.warp_data.time)/2);
                stds = -2:2;
                srvf = true;
            end

            if var_exp > 1
                error('var_exp is greater than 1')
            end
            obj.id = id;
            fn = obj.warp_data.fn;
            time = obj.warp_data.time;
            qn = obj.warp_data.qn;
            q0 = obj.warp_data.q0;
            gam = obj.warp_data.gam;

            [M, ~] = size(qn);
            
            % set up for fPCA in q-space
            mq_new = mean(qn,2);
            id = round(length(time)/2);
            if srvf
                m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));  % scaled version
                obj.mqn = [mq_new; mean(m_new)];
                qn1 = [qn; m_new];
            else
                obj.mqn = mean(fn, 2);
                qn1 = fn;
            end
            
            % calculate vector space of warping functions
            if obj.log_der
                vec = gam_to_h(gam);
                obj.mu_psi = mean(vec, 2);
            else
                [obj.mu_psi,~,~,vec] = SqrtMean(gam);
            end
            
            % joint fPCA
            f1 = @(x)find_C(x, qn1, vec, q0, no, obj.mu_psi, srvf, obj.log_der, 0.99);
            obj.C = fminbnd(f1, 0, 1e4);
            
            [qhat, gamhat, cz, Psi_q, Psi_h, sz, U, Uh, Uz] = jointfPCAd(qn1, vec, obj.C, no, obj.mu_psi, srvf, obj.log_der, var_exp);
            
            vc = C*vec;
            mv = mean(vc, 2);

            no = size(cz,2);
            % geodesic paths
            obj.stds = stds;
            obj.q_pca = zeros(M, length(stds), no);
            obj.f_pca = zeros(M, length(stds), no);
            
            for j = 1:no
                for i = 1:length(stds)
                    qhat = obj.mqn + Psi_q(:,j) * stds(i)*sqrt(sz(j));
                    vechat = Psi_h(:,j) * (stds(i)*sqrt(sz(j)))/obj.C;
                    if srvf
                        if obj.log_der
                            gamhat = h_to_gam(vechat);
                        else
                            psihat = exp_map(obj.mu_psi,vechat);
                            gamhat = cumtrapz(linspace(0,1,M), psihat.*psihat);
                        end
                        gamhat = (gamhat - min(gamhat))/(max(gamhat)-min(gamhat));
                        if (sum(vechat)==0)
                            gamhat = linspace(0,1,M);
                        end
                        
                        fhat = cumtrapzmid(time, qhat(1:M).*abs(qhat(1:M)),sign(qhat(M+1)).*(qhat(M+1)^2), id);
                        obj.f_pca(:,i,j) = warp_f_gamma(fhat, gamhat, linspace(0,1,M));
                        obj.q_pca(:,i,j) = warp_q_gamma(qhat(1:M), gamhat, linspace(0,1,M));
                    else
                        if obj.log_der
                            gamhat = h_to_gam(vechat);
                        else
                            psihat = exp_map(obj.mu_psi,vechat);
                            gamhat = cumtrapz(linspace(0,1,M), psihat.*psihat);
                        end
                        gamhat = (gamhat - min(gamhat))/(max(gamhat)-min(gamhat));
                        if (sum(vechat)==0)
                            gamhat = linspace(0,1,M);
                        end
                        
                        obj.f_pca(:,i,j) = warp_f_gamma(qhat, gamhat, linspace(0,1,M));
                        obj.q_pca(:,i,j) = f_to_srvf(obj.f_pca(:,i,j), linspace(0,1,M));
                    end
                end
            end
            
            obj.coef = cz(:, 1:no);
            obj.latent = s;
            obj.srvf = srvf;
            obj.latent = sz(1:no);
            obj.U_q = Psi_q;
            obj.U_h = Psi_h;
            obj.vec = vec;
            obj.no = no;
            obj.qn1 = qn1;
            obj.U = U;
            obj.U1 = Uh;
            obj.Uz = Uz;
            obj.mv = mv;
            
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

            if obj.srvf
                m_new = sign(fn(obj.id, :)) .* sqrt(abs(fn(obj.id, :)));
                qn1 = [qn; m_new];
            else
                qn1 = fn;
            end

            vec = zeros(M, n);
            psi = zeros(M, n);
            time = linspace(0,1,M);
            binsize = mean(diff(time));
            for i = 1:n
                if obj.log_der
                    out = gam_to_h(gam(:,i));
                else
                    psi(:, i) = sqrt(gradietn(gam(:, i), binsize));
                    [out, ~] = inv_exp_map(obj.mu_psi, psi(:,i));
                end
                vec(:, i) = out;
            end

            c = (qn1 - obj.mqn)' * obj.U;
            cv = (C*vec - obj.mv)' * obj.U1;

            Xi = [c cv];

            cz = Xi * obj.Uz;

            obj.new_coef = cz;
            obj.new_qn1 = qn1;
            obj.new_h = vec;
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
                        if i == idx
                            color = [0, 0, 0];
                        else
                            color = cl(i);
                        end
                        plot(time, obj.q_pca(:,i,k), 'Color', color, 'linewidth', 2); hold on;
                    end
                    title(['q domain: PD ' num2str(k)], 'fontsize', 14);
                    subplot(2,3,k1+3);
                    for i = 1:length(ci)
                        if i == idx
                            color = [0, 0, 0];
                        else
                            color = cl(i);
                        end
                        plot(time, obj.f_pca(:,i,k), 'Color', color, 'linewidth', 2); hold on;
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

function [qhat, gamhat, cz, Psi_q, Psi_h, sz, U, Uh, Uz] = jointfPCAd(qn, vec, C, m, mu_psi, srvf, log_der, var_exp)
[M, N] = size(qn);

% run univariate fpca
% q space
K = cov(q.');

mqn = mean(qn, 2);

[U,S,~] = svd(K);

cumm_coef = cumsum(diag(S))/sum(diag(S));
no_q = find(cumm_coef >= var_exp, 'first');
if isempty(no_q)
    no_q = 1;
end

c = (qn - mqn)' * U;
c = c(:, 1:no_q);
U = U(:, 1:no_q);

% warping functions
cvec = C*vec;
mh = mean(cvec, 2);
Kv = cov(cvec');

[Uv, Sv, ~] = svd(Kv);

cumm_coef = cumsum(diag(Sv)) / sum(diag(Sv));
no_v = find(cumm_coef >= var_exp, 'first');
if isemptyu(no_v)
    no_v = 1;
end

cv = (cvec - mh)' * Uv;
cv = cv(:, 1:no_v);
Uv = Uv(:, 1:no_v);

% run multivariate fPCA
Xi = [c cv];
Z = 1./(size(Xi,1)-1) * Xi' * Xi;

Uz, Sz, Vz = svd(Z);
cz = Xi * Uz;

Psi_q = U * Uz(1:no_q,:);
Psi_v = Uv * Uz((no_q+1):end, :);

vhat = Psi_h * cz';
if srvf
    gamhat = zeros(M-1,N);
else
    gamhat = zeros(M,N);
end
if log_der
    gamhat = h_to_gam(vhat/C);
else
    for ii = 1:N
        psihat = exp_map(mu_psi,vhat(:,ii)/C);
        if srvf
            gamhat(:,ii) = cumtrapz(linspace(0,1,M-1), psihat.*psihat);
        else
            gamhat(:,ii) = cumtrapz(linspace(0,1,M), psihat.*psihat);
        end
    end
end

qhat = Psi_q * cz' + mqn;

end

function [out] = find_C(C, qn, vec, q0, m, mu_psi, srvf, log_der, var_exp)
[qhat, gamhat, ~, ~, ~, ~, ~, ~, ~] = jointfPCAd(qn, vec, C, m, mu_psi, srvf, log_der, var_exp);
[M, N] = size(qn);

d = zeros(1,N);
for i = 1:N
    if srvf
        time = linspace(0,1,M-1);
        tmp = warp_q_gamma(qhat(1:(M-1),i), invertGamma(gamhat(:,i)), time);
    else
        time = linspace(0,1,M);
        tmp = warp_q_gamma(qhat(:,i), invertGamma(gamhat(:,i)), time);
    end
    d(i) = trapz(time,(tmp-q0(:,i)).^2);
end

out = mean(d);

end

