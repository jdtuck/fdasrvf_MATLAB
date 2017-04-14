function gamI = SqrtMeanInverse(gam)

  [n,T] = size(gam);
  time = linspace(0,1,T);

  psi = zeros(n,T);
  binsize = mean(diff(time));
  for i=1:n
      psi(i,:) = sqrt(gradient(gam(i,:),binsize));
  end

%% Find direction
mnpsi = mean(psi);
dqq = sqrt(sum((psi' - mnpsi'*ones(1,n)).^2,1));
[~, min_ind] = min(dqq);
mu = psi(min_ind,:);
maxiter = 501;
lvm = zeros(1,maxiter);
vec = zeros(n,T);
stp = .3;
iter = 1;

for i = 1:n
    vec(i,:) = inv_exp_map(mu,psi(i,:));
end
vbar=mean(vec);
lvm(iter) = L2norm(vbar);

while (lvm(iter) > 0.00000001 && iter<maxiter)
    mu = exp_map(mu, stp*vbar);
    iter = iter + 1;
    for i = 1:n
        vec(i,:) = inv_exp_map(mu,psi(i,:));
    end
    vbar=mean(vec);
    lvm(iter) = L2norm(vbar);
end

gam_mu = cumtrapz(time,mu.^2);

gam_mu = (gam_mu-min(gam_mu))/(max(gam_mu)-min(gam_mu));
gamI = invertGamma(gam_mu);
