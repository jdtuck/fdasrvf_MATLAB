function out = pcaTB(f, time, a, p, no, option)
% PCATB Computes functional tolerance bounds
% -------------------------------------------------------------------------
% Computes tolerance bounds for function data containing
% phase and amplitude variation using elasatic fPCA
%
% Usage: out = pcaTB(f, time, a, p, no)
%        out = pcaTB(f, time, a, p, no, option)
%
% Inputs:
% f (M,N): matrix defining N functions of M samples
% time : time vector of length M
% a: confidence level
% p: coverage level
% no: number of principal components
%
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting
% option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
% option.w = 0.0; % BFGS weight
% option.MaxItr = 20;  % maximum iterations
%
% Outputs:
% struct containing
% warp: alignment data from time_warping
% pca: functional pca from jointFPCA
% tol: tolerance factor

if nargin < 7
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
    option.method = 'DP1';
    option.w = 0.0;
    option.MaxItr = 20;
end

%% Align Data
out_warp = fdawarp(f,time);
out_warp = out_warp.time_warping_median(0,option);

%% Calculate pca
out_pca = jointFPCA(out_warp, no, false);

%% Calculate TB
tol = mvtol_region(out_pca.coef, a, p, 100000);

out.warp = out_warp;
out.pca = out_pca;
out.tol = tol;
