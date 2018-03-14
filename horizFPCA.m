function hfpca = horizFPCA(out_warp,no,showplot)
% HORIZFPCA Horizontal Functional Principal Component Analysis
% -------------------------------------------------------------------------
% This function calculates vertical functional principal component analysis
% on aligned data
%
% Usage: hfpca = horizFPCA(out_warp,no)
%        hfpca = horizFPCA(out_warp,no,showplot)
%
% Inputs:
% out_warp: structure from time_warping of aligned data
% no: number of principal components to extract
% showplot: show plots of principal directions (e.g, true)
%
% Outputs:
% Returns a structure containing
% gam_pca: warping functions principal directions
% psi_pca: srvf principal directions
% latent: latent values
% U: eigenvectors
% vec: shooting vectors
% mu: Karcher Mean

gam = out_warp.gam;

if nargin<3
    showplot = 1;
end

% option.showplot = 1;
[mu,~,vec] = SqrtMean(gam.');

K = cov(vec);

[U,S,~] = svd(K);
Sig = diag(S);
T = size(vec,2);

% Parameters
tau = 1:5;

% TFPCA
psi_pca = zeros(length(tau),length(mu),no);
gam_pca = zeros(length(tau),length(mu),no);
v = zeros(5,length(mu),3);
for j=1:no      % three components
    for k=tau   % -2, -1, 0, 1, 2 std from the mean
        v(k,:,j) = (k-3)*sqrt(Sig(j))*U(:,j)';
        vn = norm(v(k,:,j))/sqrt(T);
        if vn < 0.0001
            psi_pca(k,:,j) = mu;
        else
            psi_pca(k,:,j) = cos(vn)*mu.' + sin(vn)*v(k,:,j)/vn;
        end
        gam0 = cumtrapz(linspace(0,1,size(gam,1)),psi_pca(k,:,j).*psi_pca(k,:,j));
        gam_pca(k,:,j) = (gam0-gam0(1))/(gam0(end)-gam0(1));
    end
end

% coefficients
c = zeros(size(gam,2),no);
vec = vec.';
vm = mean(vec,2);
for jj = 1:no
    for ii = 1:size(gam,2)
        c(ii,jj) = (vec(:,ii)-vm).'*U(:,jj);
    end
end

hfpca.gam_pca = gam_pca;
hfpca.psi_pca = psi_pca;
hfpca.latent = Sig;
hfpca.U = U;
hfpca.coef = c;
hfpca.vec = vec;
hfpca.mu = mu;

if showplot
    cl = 'rbgmc';
    [~, ~, p1] = size(gam_pca);
    num_plot = ceil(p1/3);
    j = 1;
    for ii = 1:num_plot
        if (j > size(gam_pca,3))
            break
        end
        figure;
        for j1 = 1:3
            j = j1 + (ii-1)*3;
            if (j > size(gam_pca,3))
                break
            end
            subplot(1,3,j1);
            for k = 1:length(tau)
                plot(linspace(0,1,T), gam_pca(k,:,j), cl(k), 'linewidth', 2); hold on;
            end
            axis([0 1 0 1]);
            title(['PD ' num2str(j)], 'fontsize', 14);
        end
    end
    cumm_coef = 100*cumsum(Sig)./sum(Sig);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
end
