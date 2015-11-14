%% Original Data
clear; close all;

load toy_data;

% Align the data
[fn,qn,q0,fmean,mqn,gam,~,~] = time_warping(f,t);

C = cov(f');
mf = mean(f,2);
K = 10;  % number of samples
f_s = mvnrnd(mf', C, K)';

figure(10)
plot(t, f_s, 'LineWidth',2);
title('Random sampling on original data', 'fontsize', 14);

%% Separated and Warped Data
% sampling from the estimated model

[M, N] = size(f);
x0 = f;
x = fn;
binsize = mean(diff(t));

% compute mean and covariance in q-domain
q_new = qn;
mq_new = mqn;
m_new = sign(fn(1,:)).*sqrt(abs(fn(1,:)));  % scaled version
mqn2 = [mq_new; mean(m_new)];
C = cov([q_new;m_new]');

K = 10;  % number of samples
q_s = mvnrnd(mqn2', C, K)';

% compute the correspondence to the original function domain
for k = 1:K
    x_s(:,k) = (sign(q_s(end,k)).*(q_s(end,k).^2))...
        +cumtrapz(t,q_s(1:end-1,k).*abs(q_s(1:end-1,k)));
end

% random warping generation
rgam = randomGamma(gam.',K);


%%%% sort functions and warpings
mx = max(x_s);
[ignore, seq1] = sort(mx);

% compute the psi-function
psi = zeros(K,size(rgam,2));
for i = 1:K
    psi(i,:) = gradient(rgam(i,:), binsize)./sqrt(abs(gradient(rgam(i,:), binsize))+eps);
    ip(i) = ones(1,M)*psi(i,:)'/M;
    len(i) = acos(ones(1,M)*psi(i,:)'/M);
end
[ignore, seq2] = sort(len);

% combine x-variability and y-variability
clear x_c;
for k = 1:K
    x_c(:,k) = interp1((0:M-1)/(M-1), x_s(:,seq1(k)), invertGamma(rgam(seq2(k),:)')');
    while sum(isnan(x_c(:,k))) >= 1
        rgam2 = randomGamma(gam',1);
        x_c(:,k) = interp1((0:M-1)/(M-1), x_s(:,k), invertGamma(rgam2'));
    end
end

% plot the sampled functions
figure(11)
z = plot((0:M-1)/M,rgam, 'LineWidth',2);
% set(gca, 'xtick', 0:0.5:1, 'ytick', 0:0.2:1, 'fontsize', 20);
axis equal;
axis([0 1 0 1]);
grid minor;
title('Random sampling on x-variability', 'fontsize', 14);

figure(12)
plot(t, x_s, 'LineWidth',2);
% axis([0 1 0 9]);
% set(gca, 'xtick', 0:0.5:1, 'ytick', 4:4:8, 'fontsize', 20);
title('Random sampling on y-variability', 'fontsize', 14);

figure(13)
plot(t, x_c, 'LineWidth',2);
% axis([0 1 0 9]);
% set(gca, 'xtick', 0:0.5:1, 'ytick', 4:4:8, 'fontsize', 20);
title('Random sampling on both variabilities', 'fontsize', 14);

save model_results_toy2 x_s x_c rgam f_s