function hfpca = horizFPCA(out_warp,no,option)
% Horizontal Functional Principal Component Analysis
%
% This function calculates vertical functional principal component analysis
% on aligned data
%
% @param warp_data structure from \link{time_warping} of aligned data
% @param no number of principal components to extract
% @param showplot show plots of principal directions (default = T)
% @return Returns a structure containing
% \item{gam_pca}{warping functions principal directions}
% \item{psi_pca}{srvf principal directions}
% \item{latent}{latent values}
% \item{U}{eigenvectors}
% \item{vec}{shooting vectors}
% \item{mu}{Karcher Mean}

gam = out_warp.gam;

if nargin<3
    option.showplot = 1;
end

% option.showplot = 1;
[mu,~,vec] = SqrtMean(gam');

K = cov(vec);

[U,S,~] = svd(K);
Sig = diag(S);
T = size(vec,2) + 1;

% Parameters
tau = 1:5;

% TFPCA
psi_pca = zeros(length(tau),length(mu),no);
gam_pca = zeros(length(tau),length(mu)+1,no);
v = zeros(5,length(mu),3);
for j=1:no      % three components
    for k=tau   % -2, -1, 0, 1, 2 std from the mean
        v(k,:,j) = (k-3)*sqrt(Sig(j))*U(:,j)';
        vn = norm(v(k,:,j))/sqrt(T);
        if vn < 0.0001
            psi_pca(k,:,j) = mu;
        else
            psi_pca(k,:,j) = cos(vn)*mu + sin(vn)*v(k,:,j)/vn;
        end
        gam_pca(k,:,j) = [0 cumsum(psi_pca(k,:,j).*psi_pca(k,:,j))]/T;
    end
end

hfpca.gam_pca = gam_pca;
hfpca.psi_pca = psi_pca;
hfpca.latent = Sig;
hfpca.U = U;
hfpca.coef = c;
hfpca.vec = vec;
hfpca.mu = mu;

if option.showplot
    cl = 'rbgmc';
    [m1, n1, p1] = size(q_pca);
    num_plot = ceil(p1/3);
    cnt = 1;
    for ii = 1:num_plot
        figure;
        for j1 = 1:3
            j = j1 + (ii-1)*3
            subplot(1,3,j1);
            for k = 1:length(tau)
                plot(linspace(0,1,T), gam_pca(k,:,j), cl(j), 'linewidth', 2); hold on;
            end
            axis([0 1 0 1]);
            title(['PD ' num2str(j)], 'fontsize', 14);
        end
    end
    cumm_coef = 100*cumsum(Sig)./sum(Sig);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
end
