function gamI = SqrtMeanInverse(gam)
 
[n,T] = size(gam);
dT = 1/(T-1);
psi = zeros(n,T-1);
for i=1:n
    psi(i,:) = sqrt(diff(gam(i,:))/dT+eps);
end

%% Find direction
mnpsi = mean(psi);
dqq = sqrt(sum((psi' - mnpsi'*ones(1,n)).^2,1));
[~, min_ind] = min(dqq);
mu = psi(min_ind,:);
t = 1;
clear vec;
maxiter = 20;
lvm = zeros(1,maxiter);
vec = zeros(n,T-1);
for iter = 1:maxiter
    for i=1:n
        v = psi(i,:) - mu;
        dot1 = simps(linspace(0,1,T-1),mu.*psi(i,:));
        if dot1 > 1
            dot_limited = 1;
        elseif dot1 < -1
            dot_limited = -1;
        else
            dot_limited = dot1;
        end
        len = acos(dot_limited);
        if len > 0.0001            
            vec(i,:) = (len/sin(len))*(psi(i,:) - cos(len)*mu);
        else
            vec(i,:) = zeros(1,T-1);
        end
        
    end
    vm = mean(vec);
    lvm(iter) = sqrt(sum(vm.*vm)*dT);
    mu = cos(t*lvm(iter))*mu + (sin(t*lvm(iter))/lvm(iter))*vm;
    if lvm(iter) < 1e-6 || iter >= maxiter
        break
    end
end

gam_mu = [0 cumsum(mu.*mu)]/T;
gam_mu = (gam_mu-min(gam_mu))/(max(gam_mu)-min(gam_mu));
gamI = invertGamma(gam_mu);
