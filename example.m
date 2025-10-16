setup_paths

load data/simu_data.mat

obj = fdawarp(f, t);
obj = obj.time_warping();
obj.plot()
