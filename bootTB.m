function [amp, ph] = bootTB(f, time, a, p, B, var_exp, option)
% BOOTB Computes functional tolerance bounds
% -------------------------------------------------------------------------
% function computes tolerance bounds for function data containing
% phase and amplitude variation using bootstrap sampling
%
% Usage: [amp, ph] = bootTB(f, time, a, p, B, var_exp)
%        [amp, ph] = bootTB(f, time, a, p, B, var_exp, option)
%
% Inputs:
% f (M,N): matrix defining N functions of M samples
% time : time vector of length M
% a: confidence level (default = 0.05)
% p: coverage level (default = 0.99)
% B: number of boostrap samples (default = 1000)
% var_exp: compute number of pcs based on value percent variance explained (default = 0.99)
%
% default options
% option.parallel = 1; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines whether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting
% option.method = 'DP1'; % optimization method (DP, DP2, SIMUL, RBFGS)
% option.w = 0.0; % BFGS weight
% option.MaxItr = 20;  % maximum iterations
%
% Outputs:
% two struct containing
% amp:
%   median_y: median function
%   Q1: First quartile
%   Q3: Second quartile
%   Q1a: Lower amplitude TB based on alpha
%   Q3a: Upper amplitude TB based on alpha
%   minn: minimum extreme function
%   maxx: maximum extreme function
%   outlier_index: indexes of outlier functions
%   fmedian: median function
%   plt: surface plot
% ph:
%   median_x: median warping function
%   Q1: First quartile
%   Q3: Second quartile
%   Q1a: Lower phase TB based on alpha
%   Q3a: Upper phase TB based on alpha
%   minn: minimum extreme function
%   maxx: maximum extreme function
%   outlier_index: indexes of outlier functions
%   plt: surface plot


arguments
    f
    time
    a = 0.05;
    p = 0.99;
    B = 500;
    var_exp = 0.99;
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.method = 'DP1';
    option.spl = true;
    option.MaxItr = 20;
end

[M, ~] = size(f);

%% Align Data
out_med = fdawarp(f,time);
out_med = out_med.time_warping_median(0,'parallel',option.parallel,'closepool',option.closepool,'smooth',option.smooth, ...
    'sparam', option.sparam, 'method', option.method, 'spl', option.spl, 'MaxItr', option.MaxItr);

%% Calculate CI
if isempty(gcp('nocreate'))
    % prompt user for number threads to use
    nThreads = input('Enter number of threads to use: ');
    if nThreads > 1
        parpool(nThreads);
    elseif nThreads > 12 % check if the maximum allowable number of threads is exceeded
        while (nThreads > 12) % wait until user figures it out
            fprintf('Maximum number of threads allowed is 12\n Enter a number between 1 and 12\n');
            nThreads = input('Enter number of threads to use: ');
        end
        if nThreads > 1
            parpool(nThreads);
        end
    end
end
fprintf('Bootstrap Sampling. \n');
obj = ProgressBar(B, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd, ...
    'Title', 'Progress' ...
    );

obj.setup([], [], []);

bootlwr_amp = zeros(M,B);
bootupr_amp = zeros(M,B);
bootlwr_ph = zeros(M,B);
bootupr_ph = zeros(M,B);
parfor (k = 1:B)
    samples = joint_gauss_model(out_med, 100, var_exp);
    obja = ampbox(samples);
    amp_b = obja.construct_boxplot(1-p,.3);
    objp = phbox(samples);
    ph_b = objp.construct_boxplot(1-p,.3);
    bootlwr_amp(:,k) = amp_b.Q1a;
    bootupr_amp(:,k) = amp_b.Q3a;
    bootlwr_ph(:,k) = ph_b.Q1a;
    bootupr_ph(:,k) = ph_b.Q3a;
    updateParallel([], pwd);
end
obj.release();

if option.closepool == 1
    delete(gcp('nocreate'));
end

%% Tolerance Bounds
boot_amp = [bootlwr_amp bootupr_amp];
boot_amp_q = f_to_srvf(boot_amp, time);
boot_ph = [bootlwr_ph bootupr_ph];
boot_out = out_med;
boot_out.rsamps = 1;
boot_out.fs = boot_amp;
boot_out.qs = boot_amp_q;
boot_out.gams = boot_ph;
obj1 = ampbox(boot_out);
amp = obj1.construct_boxplot(a,.3);
obj1 = phbox(boot_out);
ph = obj1.construct_boxplot(a,.3);

