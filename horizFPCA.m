function hfpca = horizFPCA(gam,no,option)

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
psi_pca = zeros(5,length(mu),3);
gam_pca = zeros(5,length(mu)+1,3);
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

if option.showplot
    figure;
    cl = 'rbgmc';
    for j = 1:3
        subplot(1,3,j);
        for k = 1:5
            plot(linspace(0,1,T), gam_pca(k,:,j), cl(j), 'linewidth', 2); hold on;
        end
        axis([0 1 0 1]);
        title(['PD ' num2str(j)], 'fontsize', 14);
    end
    cumm_coef = 100*cumsum(Sig)./sum(Sig);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
end