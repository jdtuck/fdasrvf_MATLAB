classdef TimeWarpingConsistencyTest < matlab.unittest.TestCase
    % TIMEWARPINGCONSISTENCYTEST Aligned data must match the returned warps
    %   time_warping and time_warping_median return fn, qn and gam from the
    %   same final matching step, so fn = f o gam exactly (for smooth = 0)
    %   and, for the mean version, mqn = mean(qn).  Also checks that
    %   fdahpca.calc_fpca with the log-derivative transform gives finite,
    %   valid principal-direction warps.

    properties
        f
        t
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            S = load('data/simu_data.mat');
            testCase.f = S.f;
            testCase.t = S.t(:);
        end
    end

    methods (Test)
        function testMeanAlignmentMatchesGam(testCase)
            out = fdawarp(testCase.f, testCase.t);
            out = out.time_warping(0, MaxItr=3);
            testCase.verifyLessThan(maxWarpError(out, testCase.f), 1e-12, ...
                'fn should equal f warped by the returned gam');
            testCase.verifyLessThan(max(abs(mean(out.qn, 2) - out.mqn(:))), 1e-12, ...
                'mqn should be the mean of qn');
        end

        function testMedianAlignmentMatchesGam(testCase)
            out = fdawarp(testCase.f, testCase.t);
            out = out.time_warping_median(0, MaxItr=3);
            testCase.verifyLessThan(maxWarpError(out, testCase.f), 1e-12, ...
                'fn should equal f warped by the returned gam');
        end

        function testHorizontalPcaLogDerivative(testCase)
            out = fdawarp(testCase.f, testCase.t);
            out = out.time_warping(0, MaxItr=1);
            hpca = fdahpca(out, true);
            hpca = hpca.calc_fpca(NaN, 3);
            G = hpca.gam_pca;
            testCase.verifyTrue(all(isfinite(G(:))), ...
                'gam_pca should be finite with log_der = true');
            testCase.verifyLessThan(max(abs(G(:, 1, :) - 0), [], 'all'), 1e-12);
            testCase.verifyLessThan(max(abs(G(:, end, :) - 1), [], 'all'), 1e-12);
        end
    end
end

function e = maxWarpError(out, f)
e = 0;
for k = 1:size(f, 2)
    fw = warp_f_gamma(f(:, k), out.gam(:, k), out.time);
    e = max(e, max(abs(fw(:) - out.fn(:, k))));
end
end
