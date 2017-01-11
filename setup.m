addpath(genpath('minFunc'))
addpath(genpath('DP'))
addpath(genpath('gropt'))
addpath(genpath('mlogit_warp'))

% Dynamic Programming
fprintf('Compiling DP files...\n');
mex -outdir DP/ DP/DynamicProgrammingQ2.c DP/dp_grid.c 
mex -outdir DP/ DP/DynamicProgrammingQ.c

% minFunc (lbfgs optimizer)
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc/compiled minFunc/mex/mcholC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c

% multinomial logistic warp gradient
fprintf('Compiling mlogit files...\n');
mex -outdir mlogit_warp/ mlogit_warp/mlogit_warp.c mlogit_warp/mlogit_warp_grad.c mlogit_warp/misc_funcs.c

% compiling gropt
fprintf('Compiling gropt files...\n');
cd gropt
MyMex
cd ..
