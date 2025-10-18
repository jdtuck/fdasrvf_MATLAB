setup_paths

load data/simu_data.mat

%% Test Distance and Gamma Calc
[dy,dx] = elastic_distance(f(:,1),f(:,1),t, 0, 'DP1');
assert(dy < 1e-12, 'Amp distance not near 0')
assert(dx < 1e-7, 'Phase distance not near 0')

%% Test f warping
M = 101;
q1 = sin(linspace(0,2*pi,M));
timet = linspace(0,1,M);
gam = optimum_reparam(q1', q1', timet', 0, 'DP1'); 
q1a = warp_f_gamma(q1, gam, timet); 
assert(sum(q1(:)-q1a(:))==0,'Warping Not Identity')

%% Test q warping
M = 101;
q1 = sin(linspace(0,2*pi,M));
timet = linspace(0,1,M);
gam = optimum_reparam(q1', q1', timet', 0, 'DP1'); 
q1a = warp_q_gamma(q1, gam, timet); 
assert(sum(q1(:)-q1a(:))<=1e-15,'Warping Not Identity')

%% Test rbfgsM
M = 101;
q1 = sin(linspace(0,2*pi,M));
timet = linspace(0,1,M);
gam = optimum_reparam(q1', q1', timet', 0, 'RBFGSM'); 
q1a = warp_f_gamma(q1, gam, timet); 
assert(sum(q1(:)-q1a(:))<=1e-15,'Warping Not Identity')

%% Test Smooth
M = 101;
q1 = zeros(M,1);
q1(:,1) = sin(linspace(0,2*pi, M));
q1a = smooth_data(q1, 1);
q1b = smooth_data(q1, 1);
assert(sum(q1b(:)-q1a(:))<=1e-15,'Smoothing Not Equal')

%% Test Inv Gamma
M = 101;
gam = linspace(0, 1, M);
gami = invertGamma(gam);
assert(sum(gam(:) - gami(:)) == 0, 'invert gamma failure')

%% Test InvExpMap
M = 101;
gam = linspace(0, 1, M);
binsize = mean(diff(gam));
psi = sqrt(gradient(gam, binsize));
[out, ~] = inv_exp_map(psi, psi);
assert(sum(out(:))==0, 'InvExp Not Zero at Idenity.')

%% Test L2norm
M = 101;
gam = linspace(0, 1, M);
binsize = mean(diff(gam));
psi = sqrt(gradient(gam, binsize));
[out, theta] = inv_exp_map(psi, psi);
out1 = L2norm(out);
assert(out1==0, 'L2norm not 0.')

%% Test expmap
M = 101;
gam = linspace(0, 1, M);
binsize = mean(diff(gam));
psi = sqrt(gradient(gam, binsize));
[out, ~] = inv_exp_map(psi, psi);
out1 = exp_map(psi, out);
assert(sum(out1(:))==M, 'ExpMap Not M at Idenity.')

%% Test alignment and fpca
out = fdawarp(f, t);
out = out.time_warping(MaxItr=1);
vpca = fdavpca(out);
vpca = vpca.calc_fpca();
jpca = fdajpca(out);
jpca = jpca.calc_fpca();
hpca = fdahpca(out);
hpca = hpca.calc_fpca();
assert(length(out.time) == 101, 'Length Not Right')

%% Test amp boxplot
out = fdawarp(f, t);
out = out.time_warping_median(MaxItr=1);
amp_data = ampbox(out);
amp_data = amp_data.construct_boxplot(.05, 1);
assert(length(amp_data.Q1) == 101)

%% Test ph boxplot
out = fdawarp(f, t);
out = out.time_warping_median(MaxItr=1);
ph_data = phbox(out);
ph_data = ph_data.construct_boxplot(.05, 1);
assert(length(ph_data.Q1) == 101)

%% Test Curve Align
load data/MPEG7.mat
curve = Xdata{1,2};
[n, M] = size(curve);
K = size(Xdata,2);

beta = zeros(n, M, K);
for i = 1:K
    beta(:, :, i) = Xdata{1, i};
end

obj = fdacurve(beta, false, M);
obj = obj.karcher_mean();
obj = obj.karcher_cov();
obj = obj.shape_pca();
