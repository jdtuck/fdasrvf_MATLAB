classdef fdaspc
    %fdaspc A class to provide elastic functional statistical process control
    % -------------------------------------------------------------------------
    % This class provides PCA-based process monitoring (Hotelling T^2 and
    % SPE/Q charts) of functional data in the SRVF framework, with separate
    % models for amplitude (vertical fPCA) and phase (horizontal fPCA). The
    % monitoring statistics and control limits follow Z. Ge and Z. Song,
    % "Multivariate Statistical Process Control", Springer 2013, Sec. 2.1,
    % Eqs. (2.3)-(2.8), applied to the elastic fPCA features.
    %
    % Reference:
    %   J. D. Tucker, "Functional Statistical Process Control using Elastic
    %   Methods," Proc. JSM, August 2016.
    %   https://research.tetonedge.net/Papers/JSM16.pdf
    %
    % Usage:  obj = fdaspc(warp_data)
    %         obj = obj.train()
    %         [obj, out] = obj.monitor(f_new)
    %
    % where:
    %   warp_data - fdawarp object of aligned normal-operating-condition data
    %   f_new     - (M,N) matrix of N new functions sampled at warp_data.time
    %
    %
    % fdaspc Properties:
    %   warp_data - fdawarp class with alignment data
    %   vpca - fdavpca object (amplitude model)
    %   hpca - fdahpca object (phase model)
    %   alpha - significance level of the control limits
    %   t2limit - T^2 limit formula, 'standard' or 'book'
    %   amp - amplitude model: k, n, T2lim, SPElim, g, h, training T2/SPE
    %   ph - phase model: k, n, T2lim, SPElim, g, h, training T2/SPE
    %   mon - results of the last call to monitor
    %
    %
    % fdaspc Methods:
    %   fdaspc - class constructor
    %   train - fit the amplitude and phase fPCA models and control limits
    %   monitor - align new functions and compute T^2/SPE and alarms
    %   plot - control charts of the training and monitored data
    %   plot_contrib - T^2 and SPE contributions over time for one function
    %
    %
    % Author :  J. D. Tucker (JDT) <jdtuck AT sandia.gov>
    % Date   :  20-Nov-2026

    properties
        warp_data % fdawarp class with alignment data
        vpca      % fdavpca object (amplitude model)
        hpca      % fdahpca object (phase model)
        alpha     % significance level of the control limits
        t2limit   % T^2 limit formula, 'standard' or 'book'
        amp       % amplitude model and training statistics
        ph        % phase model and training statistics
        mon       % results of the last call to monitor
    end

    methods
        function obj = fdaspc(fdawarp)
            %fdaspc Construct an instance of this class
            % Input:
            %   fdawarp: fdawarp class aligned with time_warping or
            %            time_warping_median
            if (isempty(fdawarp.fn))
                error('Please align fdawarp class using time_warping!');
            end
            obj.warp_data = fdawarp;
        end

        function obj = train(obj, option)
            % TRAIN Fit the amplitude and phase monitoring models
            % -------------------------------------------------------------------------
            % Runs vertical and horizontal fPCA on the aligned data and
            % computes the Hotelling T^2 limit (Ge & Song Eq. 2.7) and SPE
            % limit (Eq. 2.8, Box chi-square approximation) for each.
            %
            % Usage: obj = obj.train()
            %        obj = obj.train(var_exp=0.9, alpha=0.05)
            %        obj = obj.train(no=[3 2])
            %
            % Inputs (name-value):
            % var_exp: cumulative variance used to choose the number of
            %          components (default 0.95); scalar, or [amp phase]
            % no: number of components; scalar or [amp phase]. When given,
            %     it overrides var_exp (default NaN)
            % alpha: significance level of the limits (default 0.01)
            % t2limit: 'standard', k(n^2-1)/(n(n-k)) F, or 'book',
            %          k(n-1)/(n-k) F (default 'standard')
            % id: point used for f(id) in the amplitude model
            %     (default midpoint)
            % log_der: use the log-derivative transform for the phase
            %          model instead of the shooting vectors (default false)
            % realign: compute the training T^2/SPE and the limits from
            %          the training functions re-aligned to the template
            %          with align_new, the same path monitor uses (default
            %          true). Re-alignment perturbs the srvfs, so limits
            %          fit to the in-sample residuals are too tight for
            %          monitored data; false uses the in-sample features
            % spl: srvf derivative used for the re-alignment (default true)
            %
            % Output:
            % fdaspc object
            arguments
                obj
                option.var_exp = 0.95;
                option.no = NaN;
                option.alpha = 0.01;
                option.t2limit {mustBeMember(option.t2limit, ["standard","book"])} = "standard";
                option.id = round(length(obj.warp_data.time)/2);
                option.log_der = false;
                option.realign = true;
                option.spl = true;
            end

            if ~(option.alpha > 0 && option.alpha < 1)
                error('alpha must be in (0, 1)');
            end
            var_exp = pair(option.var_exp);
            no = pair(option.no);
            var_exp(~isnan(no)) = NaN;

            obj.alpha = option.alpha;
            obj.t2limit = char(option.t2limit);

            % amplitude: vertical fPCA on [qn; sign(f(id))sqrt|f(id)|]
            obj.vpca = fdavpca(obj.warp_data);
            obj.vpca = obj.vpca.calc_fpca(var_exp(1), no(1), option.id);

            % phase: horizontal fPCA on the shooting vectors
            obj.hpca = fdahpca(obj.warp_data, option.log_der);
            obj.hpca = obj.hpca.calc_fpca(var_exp(2), no(2));

            % features the limits are fit to
            if option.realign
                [fn, qn, gam] = obj.warp_data.align_new(obj.warp_data.f, ...
                    spl=option.spl);
                X = amp_features(obj.vpca, fn, qn);
                V = ph_features(obj.hpca, gam);
            else
                X = amp_features(obj.vpca, obj.warp_data.fn, obj.warp_data.qn);
                V = obj.hpca.vec;
            end
            obj.amp = fit_limits(X, obj.vpca.mqn, obj.vpca.U, ...
                obj.vpca.latent, obj.alpha, obj.t2limit);
            obj.ph = fit_limits(V, obj.hpca.vm, obj.hpca.U, ...
                obj.hpca.latent, obj.alpha, obj.t2limit);

            obj.mon = [];
        end

        function [obj, out] = monitor(obj, f_new, option)
            % MONITOR Monitor new functions against the trained models
            % -------------------------------------------------------------------------
            % Aligns the new functions to the training template (same
            % lambda, method and penalty) and computes the scores, T^2, SPE,
            % alarms and per-variable contributions for amplitude and phase.
            % A fault is detected when T^2 > T2lim or SPE > SPElim.
            %
            % Usage: [obj, out] = obj.monitor(f_new)
            %
            % Inputs:
            % f_new: (M,N) matrix of N new functions sampled at
            %        warp_data.time
            % spl: srvf derivative used for the alignment fit
            %      (default true, see fdawarp.align_new)
            %
            % Outputs:
            % obj: fdaspc object with mon filled in
            % out: struct with fields
            %   fn, qn, gam - aligned functions, srvfs and warps
            %   amp, ph - structs with T (scores), T2, SPE, alarmT2,
            %             alarmSPE, alarm, contT2, contSPE (rows are new
            %             functions; contributions are by time point, the
            %             amplitude ones have the f(id) term last)
            %   alarm - (N,1) logical, any amplitude or phase alarm
            arguments
                obj
                f_new
                option.spl = true;
            end
            if isempty(obj.amp)
                error('Please train the fdaspc object first!');
            end

            [fn, qn, gam] = obj.warp_data.align_new(f_new, spl=option.spl);

            % amplitude
            X = amp_features(obj.vpca, fn, qn);
            amp_out = monitor_stats(X, obj.vpca.mqn, obj.vpca.U, ...
                obj.vpca.latent, obj.amp);

            % phase
            V = ph_features(obj.hpca, gam);
            ph_out = monitor_stats(V, obj.hpca.vm, obj.hpca.U, ...
                obj.hpca.latent, obj.ph);

            out = struct('fn', fn, 'qn', qn, 'gam', gam, ...
                'amp', amp_out, 'ph', ph_out, ...
                'alarm', amp_out.alarm | ph_out.alarm);
            obj.mon = out;
        end

        function plot(obj)
            % PLOT control charts for amplitude and phase
            % -------------------------------------------------------------------------
            % T^2 and SPE charts of the training data, followed by the data
            % from the last call to monitor; dashed lines are the control
            % limits and red points are alarms.
            %
            % Usage: obj.plot()
            if isempty(obj.amp)
                error('Please train the fdaspc object first!');
            end
            mods = {obj.amp, obj.ph};
            names = {'Amplitude', 'Phase'};
            stats = {'T2', 'SPE'};
            labels = {'T^2', 'SPE'};
            if isempty(obj.mon)
                news = {[], []};
            else
                news = {obj.mon.amp, obj.mon.ph};
            end

            figure;
            tiledlayout(2, 2);
            for a = 1:2
                for s = 1:2
                    nexttile;
                    tr = mods{a}.(stats{s});
                    lim = mods{a}.([stats{s} 'lim']);
                    ntr = length(tr);
                    plot(1:ntr, tr, '.', 'Color', [0.5 0.5 0.5]); hold on;
                    if ~isempty(news{a})
                        nw = news{a}.(stats{s});
                        idx = ntr + (1:length(nw));
                        al = nw > lim;
                        plot(idx, nw, 'b.-');
                        plot(idx(al), nw(al), 'ro', 'MarkerFaceColor', 'r');
                        xline(ntr + 0.5, ':');
                    end
                    yline(lim, 'k--', 'LineWidth', 1.5);
                    hold off;
                    xlabel('Observation');
                    ylabel(labels{s});
                    title(sprintf('%s %s (k = %d)', names{a}, labels{s}, mods{a}.k));
                end
            end
        end

        function plot_contrib(obj, idx)
            % PLOT_CONTRIB contributions over time for a monitored function
            % -------------------------------------------------------------------------
            % Plots the T^2 and SPE contributions of each time point for
            % the amplitude and phase models, to locate where in time a
            % fault occurs.
            %
            % Usage: obj.plot_contrib(idx)
            %
            % Input:
            % idx: index of the function passed to the last monitor call
            if isempty(obj.mon)
                error('Please call monitor first!');
            end
            t = obj.warp_data.time;
            figure;
            tiledlayout(2, 2);
            nexttile; plot(t, obj.mon.amp.contT2(idx, 1:end-1), 'LineWidth', 1.5);
            title('Amplitude T^2 contribution'); xlabel('t');
            nexttile; plot(t, obj.mon.amp.contSPE(idx, 1:end-1), 'LineWidth', 1.5);
            title('Amplitude SPE contribution'); xlabel('t');
            nexttile; plot(t, obj.mon.ph.contT2(idx, :), 'LineWidth', 1.5);
            title('Phase T^2 contribution'); xlabel('t');
            nexttile; plot(t, obj.mon.ph.contSPE(idx, :), 'LineWidth', 1.5);
            title('Phase SPE contribution'); xlabel('t');
        end
    end
end

% ======================================================================
function v = pair(x)
% scalar -> [x x], keep a 2-vector
if isscalar(x)
    v = [x x];
elseif numel(x) == 2
    v = x(:).';
else
    error('expected a scalar or a 2-element [amp phase] vector');
end
end

function X = amp_features(vpca, fn, qn)
% columns [qn; sign(f(id))sqrt|f(id)|] as in fdavpca
m = sign(fn(vpca.id,:)) .* sqrt(abs(fn(vpca.id,:)));
X = [qn; m];
end

function V = ph_features(hpca, gam)
% shooting vectors at the Karcher mean (same map as SqrtMean /
% fdahpca.project), or the log-derivative transform
if hpca.log_der
    V = gam_to_h(gam);
else
    [M, n] = size(gam);
    binsize = 1/(M-1);
    V = zeros(M, n);
    for i = 1:n
        psi = sqrt(max(gradient(gam(:,i), binsize), 0));
        V(:,i) = inv_exp_map(hpca.mu, psi);
    end
end
end

function S = scores(X, mu, U, latent)
% Ge & Song Eqs. (2.3)-(2.6) and contributions; X is (p,N)
Xc = (X - mu(:)).';                  % N x p
T = Xc * U;                          % scores
E = Xc - T * U.';                    % residuals
lam = latent(:).';
T2 = sum(T.^2 ./ lam, 2);
SPE = sum(E.^2, 2);
D = U ./ sqrt(lam);
contT2 = (Xc * D) * D.' .* Xc;       % rows sum to T2
contSPE = E.^2;                      % rows sum to SPE
S = struct('T', T, 'T2', T2, 'SPE', SPE, ...
    'contT2', contT2, 'contSPE', contSPE);
end

function mdl = fit_limits(X, mu, U, latent, alpha, t2limit)
% control limits from the training data, Ge & Song Eqs. (2.7)-(2.8)
n = size(X, 2);
k = size(U, 2);
if k >= n
    error('fdaspc:kTooLarge', ...
        'Need fewer components than functions for the T^2 limit (k = %d, n = %d).', k, n);
end
S = scores(X, mu, U, latent);

Fcrit = finv(1 - alpha, k, n - k);
switch t2limit
    case 'book'
        T2lim = k * (n - 1) / (n - k) * Fcrit;
    otherwise
        T2lim = k * (n^2 - 1) / (n * (n - k)) * Fcrit;
end

mSPE = mean(S.SPE);
vSPE = var(S.SPE);
if vSPE < eps
    g = eps; h = 1;
else
    g = vSPE / (2 * mSPE);
    h = 2 * mSPE^2 / vSPE;
end
SPElim = g * chi2inv(1 - alpha, h);

mdl = struct('k', k, 'n', n, 'T2lim', T2lim, 'SPElim', SPElim, ...
    'g', g, 'h', h, 'T2', S.T2, 'SPE', S.SPE);
end

function out = monitor_stats(X, mu, U, latent, mdl)
S = scores(X, mu, U, latent);
out = S;
out.alarmT2 = S.T2 > mdl.T2lim;
out.alarmSPE = S.SPE > mdl.SPElim;
out.alarm = out.alarmT2 | out.alarmSPE;
out.T2lim = mdl.T2lim;
out.SPElim = mdl.SPElim;
end
