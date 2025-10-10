% setup environment
setup_paths

% Dynamic Programming
fprintf('Compiling DP files...\n');
mex -outdir DP/ DP/DynamicProgrammingQ2.c DP/dp_grid.c DP/dp_nbhd.c
mex -outdir DP/ DP/DynamicProgrammingQ.c

% Bayesian
fprintf('Compiling Bayesian files...\n');
if ismac
    mex -ld_classic -outdir armadillo_cpp/ armadillo_cpp/bcalcY.cpp
    mex -ld_classic -outdir armadillo_cpp/ armadillo_cpp/bcuL2norm2.cpp
    mex -ld_classic -outdir armadillo_cpp/ armadillo_cpp/trapzCpp.cpp
    mex -ld_classic -outdir armadillo_cpp/ armadillo_cpp/border_l2norm.cpp
else
    mex -outdir armadillo_cpp/ armadillo_cpp/bcalcY.cpp
    mex -outdir armadillo_cpp/ armadillo_cpp/bcuL2norm2.cpp
    mex -outdir armadillo_cpp/ armadillo_cpp/trapzCpp.cpp
    mex -outdir armadillo_cpp/ armadillo_cpp/border_l2norm.cpp
end

% minFunc (lbfgs optimizer)
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc/compiled minFunc/mex/mcholC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c

% multinomial logistic warp gradient
fprintf('Compiling mlogit files...\n');
mex -outdir mlogit_warp/ mlogit_warp/mlogit_warp.c mlogit_warp/mlogit_warp_grad.c mlogit_warp/misc_funcs.c

% Bayesian
fprintf('Compiling rbfgs files...\n');
if ispc
    libDir = fullfile(matlabroot,"/extern/lib/win64/microsoft");    %Path to LAPAC and BLAS libs
    mexCmd = sprintf('mex LINKFLAGS=''$LINKFLAGS -ld_classic'' -outdir armadillo_cpp/ armadillo_cpp/c_rlbfgs.cpp -L"%s" -lmwlapack -lmwblas', libDir);
    eval(mexCmd);
elseif ismac
    mex -ld_classic -outdir armadillo_cpp/ armadillo_cpp/c_rlbfgs.cpp -llapack -lblas
else
    mex -outdir armadillo_cpp/ armadillo_cpp/c_rlbfgs.cpp -llapack -lblas
end
