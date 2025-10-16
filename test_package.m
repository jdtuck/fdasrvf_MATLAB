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
