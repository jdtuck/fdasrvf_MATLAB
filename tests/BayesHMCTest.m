classdef BayesHMCTest < matlab.unittest.TestCase
    % BAYESHMCTEST Bayesian HMC alignment (pairwise and group)
    %   Checks the exact likelihood gradient used by the infinite-dimensional
    %   HMC sampler against finite differences for both warp representations,
    %   that the compiled kernel (when built) matches the MATLAB code, and
    %   that short pairwise and group runs return valid warps close to the
    %   dynamic programming solution.

    properties
        f
        t
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            S = load('simu_data.mat');
            testCase.f = S.f;
            testCase.t = S.t(:);
        end
    end

    methods (Test)
        function testGradientMatchesFiniteDifference(testCase)
            [Phi, q1, q2] = smallProblem();
            rng(1);
            c = 0.2 * randn(size(Phi, 2), 1);
            for rep = 0:1
                [~, g] = bayes_warp_geom(c, Phi, rep, q1, q2, 0.1);
                fd = zeros(size(c));
                e = 1e-7;
                for i = 1:numel(c)
                    cp = c; cp(i) = cp(i) + e;
                    cm = c; cm(i) = cm(i) - e;
                    fd(i) = (bayes_warp_geom(cp, Phi, rep, q1, q2, 0.1) ...
                        - bayes_warp_geom(cm, Phi, rep, q1, q2, 0.1)) / (2*e);
                end
                testCase.verifyLessThan(max(abs(fd - g)) / max(abs(fd)), 1e-5);
            end
        end

        function testShootingOutsideInjectivityRadius(testCase)
            [Phi, q1, q2] = smallProblem();
            c = zeros(size(Phi, 2), 1);
            c(1) = 4; % |v| = 4 > pi
            testCase.verifyEqual(bayes_warp_geom(c, Phi, 0, q1, q2, 0.1), Inf);
        end

        function testCoefGamRoundTrip(testCase)
            [Phi, ~, ~] = smallProblem();
            rng(2);
            c = 0.1 * randn(size(Phi, 2), 3);
            for rep = 0:1
                gam = bayes_coef_to_gam(c, Phi, rep);
                testCase.verifyEqual(gam(1, :), zeros(1, 3));
                testCase.verifyEqual(gam(end, :), ones(1, 3));
                testCase.verifyTrue(all(diff(gam) >= 0, 'all'));
                c2 = bayes_gam_to_coef(gam, Phi, rep);
                testCase.verifyLessThan(max(abs(c2 - c), [], 'all'), 1e-2);
            end
        end

        function testMexMatchesMatlab(testCase)
            testCase.assumeEqual(exist('bayes_warp_geom_c', 'file'), 3, ...
                'bayes_warp_geom_c is not built');
            [Phi, q1, q2] = smallProblem();
            rng(3);
            c = 0.3 * randn(size(Phi, 2), 1);
            for rep = 0:1
                [n1, g1, s1, gam1, w1] = bayes_warp_geom_c(c, Phi, rep, q1, q2, 0.1);
                [n2, g2, s2, gam2, w2] = matlabGeom(c, Phi, rep, q1, q2, 0.1);
                testCase.verifyEqual(n1, n2, 'AbsTol', 1e-9);
                testCase.verifyEqual(s1, s2, 'AbsTol', 1e-9);
                testCase.verifyEqual(g1, g2, 'RelTol', 1e-9);
                testCase.verifyEqual(gam1, gam2, 'AbsTol', 1e-12);
                testCase.verifyEqual(w1, w2, 'AbsTol', 1e-12);
            end
        end

        function testBatchedGeomMatchesSingle(testCase)
            [Phi, q1, ~] = smallProblem();
            rng(6);
            tt = linspace(0, 1, size(Phi, 1)).';
            Q = f_to_srvf(testCase.f(:, 1:5), tt);
            c = 0.2 * randn(size(Phi, 2), 5);
            for rep = 0:1
                [n, g, s, gam, w] = bayes_warp_geom(c, Phi, rep, q1, Q, 0.1);
                for i = 1:5
                    [n1, g1, s1, gam1, w1] = bayes_warp_geom(c(:,i), Phi, rep, q1, Q(:,i), 0.1);
                    testCase.verifyEqual([n(i), s(i)], [n1, s1], 'AbsTol', 1e-10);
                    testCase.verifyEqual(g(:,i), g1, 'AbsTol', 1e-10);
                    testCase.verifyEqual(gam(:,i), gam1, 'AbsTol', 1e-14);
                    testCase.verifyEqual(w(:,i), w1, 'AbsTol', 1e-14);
                end
            end
        end

        function testBatchedHmcSamplesGaussian(testCase)
            % nll = 0: the target is the prior N(0, diag(C)), which the
            % rotation integrates exactly (every proposal is accepted)
            rng(8);
            C = [1; 0.25; 0.04];
            geom = @zeroGeom;
            c = zeros(3, 2000);
            [nll, g] = geom(c);
            for k = 1:20
                [c, nll, g, ~, acc] = infhmc_step_batch(c, nll, g, [], geom, C, 0.7, 5);
            end
            testCase.verifyTrue(all(acc));
            testCase.verifyEqual(var(c, 0, 2), C, 'RelTol', 0.1);
        end

        function testPairwiseAlignment(testCase)
            rng(4);
            y1 = testCase.f(:,1).' + 0.01 * randn(1, numel(testCase.t));
            y2 = testCase.f(:,21).' + 0.01 * randn(1, numel(testCase.t));
            [~, gdp] = pairwise_align(testCase.f(:,1), testCase.f(:,21), testCase.t);
            for rep = ["shooting", "clr"]
                evalc('out = pairwise_align_bayes_infHMCbasis(y1, y2, testCase.t.'', iter=1000, thin=1, L=10, warp_rep=rep);');
                testCase.verifyGreaterThan(mean(out.theta_accept(500:end)), 0.3);
                testCase.verifyLessThan(max(abs(out.gamma(:) - gdp(:))), 0.08);
                testCase.verifyEqual(size(out.gamma_mat), [numel(testCase.t), 501]);
            end
        end

        function testGroupAlignment(testCase)
            rng(5);
            d = fdawarp(testCase.f, testCase.t);
            evalc('d = d.time_warping(0, MaxItr=10);');
            obj = fdawarp(testCase.f, testCase.t);
            evalc('obj = obj.time_warping_bayes_infHMC(iter=600);');
            M = numel(testCase.t);
            N = size(testCase.f, 2);
            testCase.verifyEqual(size(obj.gam), [M, N]);
            testCase.verifyTrue(all(diff(obj.gam) >= -1e-12, 'all'));
            testCase.verifyLessThan(mean(abs(obj.gam - d.gam), 'all'), 0.02);
            testCase.verifyGreaterThan(mean(obj.mcmc.accept(:, 300:end), 'all'), 0.4);
            fn = zeros(M, N);
            for k = 1:N
                fn(:,k) = warp_f_gamma(testCase.f(:,k), obj.gam(:,k), testCase.t);
            end
            testCase.verifyLessThan(max(abs(fn - obj.fn), [], 'all'), 1e-12);
            testCase.verifyEqual(size(obj.mcmc.gamma_stats), [M, 2, N]);
        end
    end
end

function [nll, g, aux] = zeroGeom(c)
nll = zeros(1, size(c, 2));
g = zeros(size(c));
aux = [];
end

function [Phi, q1, q2] = smallProblem()
T = 101;
t = linspace(0, 1, T).';
Phi = bayes_warp_basis(t, 10, 'fourier');
q1 = f_to_srvf(exp(-(t - 0.4).^2 / 0.01), t);
q2 = f_to_srvf(exp(-(t - 0.6).^2 / 0.01), t);
end

function [nll, grad, SSE, gam, q2w] = matlabGeom(c, Phi, rep, q1, q2, var1)
% call the MATLAB implementation of bayes_warp_geom, bypassing the MEX
src = fileread(which('bayes_warp_geom'));
src = strrep(src, "has_mex = exist('bayes_warp_geom_c', 'file') == 3;", 'has_mex = false;');
src = strrep(src, 'function [nll, grad, SSE, gam, q2w] = bayes_warp_geom(', ...
    'function [nll, grad, SSE, gam, q2w] = bwg_matlab_ref(');
d = tempname;
mkdir(d);
fid = fopen(fullfile(d, 'bwg_matlab_ref.m'), 'w');
fwrite(fid, src);
fclose(fid);
addpath(d);
cleanup = onCleanup(@() removeDir(d));
[nll, grad, SSE, gam, q2w] = bwg_matlab_ref(c, Phi, rep, q1, q2, var1);
end

function removeDir(d)
rmpath(d);
rmdir(d, 's');
end
