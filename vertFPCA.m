function vfpca = vertFPCA(fn,t,qn,no,option)
if nargin<5
    option.showplot = 1;
end

% Parameters
coef = -2:2;
NP = 1:no;  % number of principal components
Nstd = length(coef);

% FPCA
mq_new = mean(qn,2);
m_new = sign(fn(round(length(t)/2),:)).*sqrt(abs(fn(round(length(t)/2),:)));  % scaled version
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
        f_pca(:,i,k) = cumtrapzmid(t,q_pca(1:end-1,i,k).*abs(q_pca(1:end-1,i,k)),sign(q_pca(end,i,k)).*(q_pca(end,i,k).^2));
    end
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
vfpca.c = c;
vfpca.U = U;

if option.showplot
    cl = 'rbgmc';
    figure
    for k = NP
        subplot(2,3,k);
        for i = 1:length(coef)
            plot(t, q_pca(1:end-1,i,k), cl(i), 'linewidth', 2); hold on;
        end
        title(['q domain: PD ' num2str(k)], 'fontsize', 14);
        subplot(2,3,k+3);
        for i = 1:length(coef)
            plot(t, f_pca(:,i,k), cl(i), 'linewidth', 2); hold on;
        end
        title(['f domain: PD ' num2str(k)], 'fontsize', 14);
    end
    cumm_coef = 100*cumsum(s)./sum(s);
    figure
    plot(cumm_coef);title('Coefficient Cumulative Percentage');ylabel('Percentage');xlabel('Index')
end
