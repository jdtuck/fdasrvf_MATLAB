function rgam = randomGamma(gam,num)

[mu,psi,vec] = SqrtMean(gam);

K = cov(vec);
[U,S,V] = svd(K);
Sig = diag(S);
n = 5;
T = length(vec) + 1;

% Random Sampling
% figure(12); clf; 
% axes('FontSize',20);
% hold on;

vm = mean(vec);
for k=1:num
    
    a = randn(1,n);
    v = zeros(size(vm));
    for i=1:n
        v = v + a(i)*sqrt(Sig(i))*U(:,i)';
    end
    vn = norm(v)/sqrt(T);
    psi = cos(vn)*mu + sin(vn)*v/vn;

    gam0 = [0 cumsum(psi.*psi)]/T;
    rgam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale

    % plot([1:T]/T,rgam(k,:),'LineWidth',2);
end
% axis equal;
% axis([0 1 0 1]);
% grid minor;

