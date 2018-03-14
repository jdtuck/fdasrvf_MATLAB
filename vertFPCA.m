function vfpca = vertFPCA(out_warp,no,id,showplot)
% VERTFPCA Vertical Functional Principal Component Analysis
% -------------------------------------------------------------------------
% This function calculates vertical functional principal component analysis
% on aligned data
%
% Usage: vfpca = vertFPCA(out_warp,no)
%        vfpca = vertFPCA(out_warp,no,id)
%        vfpca = vertFPCA(out_warp,no,id,showplot)
%
% Inputs:
% warp_data: struct from time_warping of aligned data
% no: number of principal components to extract
% id: point to use for f(0) (default = midpoint)
% showplot: show plots of principal directions (default = true)
%
% Output
% Returns a struct containing
% q_pca: srvf principal directions
% f_pca: f principal directions
% latent: latent values
% coef: coefficients
% U: eigenvectors
% id: point used for f(0)
fn = out_warp.fn;
t = out_warp.time;
qn = out_warp.qn;

if nargin < 3
    id = round(length(t)/2);
    showplot = 1;
elseif nargin < 4
    showplot = 1;
end

% Parameters
coef = -2:2;
NP = 1:no;  % number of principal components
Nstd = length(coef);

% FPCA
mq_new = mean(qn,2);
m_new = sign(fn(id,:)).*sqrt(abs(fn(id,:)));  % scaled version
mqn = [mq_new; mean(m_new)];
K = cov([qn;m_new]');

[U,S,~] = svd(K);
s = diag(S);
stdS = sqrt(s);

% compute the PCA in the q domain
q_pca = zeros(length(mq_new)+1,Nstd,no);
for k = NP
    for i = 1:Nstd
        q_pca(:,i,k) = [mq_new; mean(m_new)] + coef(i)*stdS(k)*U(:,k);
    end
end

% compute the correspondence to the original function domain
f_pca = zeros(length(mq_new),Nstd,no);
for k = NP
    for i = 1:Nstd
        f_pca(:,i,k) = cumtrapzmid(t,q_pca(1:end-1,i,k).*abs(q_pca(1:end-1,i,k)),sign(q_pca(end,i,k)).*(q_pca(end,i,k).^2), id);
    end
    fbar = mean(fn,2);
    fsbar = mean(f_pca(:,:,k),2);
    err = repmat(fbar-fsbar,1,Nstd);
    f_pca(:,:,k) = f_pca(:,:,k) + err;
end

% coefficients
c = zeros(size(qn,2),no);
for jj = NP
    for ii = 1:size(qn,2)
        c(ii,jj) = dot([qn(:,ii);m_new(ii)]-mqn,U(:,jj));
    end
end

vfpca.q_pca = q_pca;
vfpca.f_pca = f_pca;
vfpca.latent = s;
vfpca.coef = c;
vfpca.id = id;
vfpca.mqn = mqn;
vfpca.time = t;
vfpca.c = c;
vfpca.U = U;

if showplot
    cl = 'rbgmc';
    [~, ~, p1] = size(q_pca);
    num_plot = ceil(p1/3);
    for ii = 1:num_plot
        if (k > size(q_pca,3))
            break
        end
        figure;
        for k1 = 1:3
            k = k1+(ii-1)*3;
            subplot(2,3,k1);
            if (k > size(q_pca,3))
                break
            end
            for i = 1:length(coef)
                plot(t, q_pca(1:end-1,i,k), cl(i), 'linewidth', 2); hold on;
            end
            title(['q domain: PD ' num2str(k)], 'fontsize', 14);
            subplot(2,3,k1+3);
            for i = 1:length(coef)
                plot(t, f_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
            end
            title(['f domain: PD ' num2str(k)], 'fontsize', 14);
        end
    end
    cumm_coef = 100*cumsum(s)./sum(s);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
end
