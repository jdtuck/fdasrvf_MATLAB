classdef DynamicProgrammingTest < matlab.unittest.TestCase
    % DYNAMICPROGRAMMINGTEST Accuracy of the DP alignments on a known warp
    %   f2 = f1 o gam_true, so aligning q2 to q1 should recover the inverse
    %   of gam_true.  The DP cost resamples q1 and q2 three times finer
    %   (scl = 3) and interpolates q2 between those points; rounding q2 to
    %   the nearest grid point at scl = 1 left the warp up to 0.63 grid cells
    %   off and the fit worse than the true warp's.

    properties
        t
        q1
        q2
        gamInv
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            M = 101;
            testCase.t = linspace(0, 1, M)';
            f = @(s) sin(3*pi*s) + 0.5*cos(7*pi*s) + 2*s;
            gamTrue = testCase.t + 0.5*sin(pi*testCase.t)/pi;
            testCase.q1 = f_to_srvf(f(testCase.t), testCase.t);
            testCase.q2 = f_to_srvf(f(gamTrue), testCase.t);
            testCase.gamInv = interp1(gamTrue, testCase.t, testCase.t);
        end
    end

    methods (Test)
        function testRecoversKnownWarp(testCase)
            gam = optimum_reparam(testCase.q1, testCase.q2, testCase.t, 0, 'DP1');
            M = length(testCase.t);
            testCase.verifyLessThan(max(abs(gam(:) - testCase.gamInv)) * (M-1), 0.5, ...
                'DP1 warp should be within half a grid cell of the true warp');
        end

        function testFitsAsWellAsTrueWarp(testCase)
            gam = optimum_reparam(testCase.q1, testCase.q2, testCase.t, 0, 'DP1');
            dy = alignError(testCase, gam(:));
            dyTrue = alignError(testCase, testCase.gamInv);
            testCase.verifyLessThan(dy, 1.1 * dyTrue, ...
                'DP1 alignment should fit about as well as the true warp');
        end

        function testDefaultScale(testCase)
            % scl defaults to 3 and can be set as a fifth argument
            q1 = testCase.q1';
            q2 = testCase.q2';
            g = DynamicProgrammingQ(q2, q1, 0, 0);
            testCase.verifyEqual(DynamicProgrammingQ(q2, q1, 0, 0, 3), g, ...
                'DynamicProgrammingQ should default to scl = 3');
            g1 = DynamicProgrammingQ(q2, q1, 0, 0, 1);
            testCase.verifyEqual(size(g1), size(g));
            testCase.verifyTrue(all(diff(g1) >= 0), 'scl = 1 warp should be monotone');
        end

        function testDPFitsAsWellAsTrueWarp(testCase)
            % "DP" resamples the SRVFs scl = 3 times finer by default; with
            % scl = 1 its fit was 1.5x worse than the true warp's
            gam = optimum_reparam(testCase.q1, testCase.q2, testCase.t, 0, 'DP');
            dy = alignError(testCase, gam(:));
            dyTrue = alignError(testCase, testCase.gamInv);
            testCase.verifyLessThan(dy, 1.1 * dyTrue, ...
                'DP alignment should fit about as well as the true warp');
        end

        function testSclPassedThrough(testCase)
            % scl = 1 gives the solvers' unresampled results
            t = testCase.t;
            q1n = testCase.q1 / norm(testCase.q1);
            q2n = testCase.q2 / norm(testCase.q2);
            gam = optimum_reparam(testCase.q1, testCase.q2, t, 0, 'DP1', 0, 0, 7, 'roughness', 1);
            g = DynamicProgrammingQ(q2n', q1n', 0, 1, 1);
            testCase.verifyEqual(gam(:), (g(:) - g(1)) / (g(end) - g(1)), 'AbsTol', 1e-14);
            gam = optimum_reparam(testCase.q1, testCase.q2, t, 0, 'DP', 0, 0, 7, 'roughness', 1);
            [G, T] = DynamicProgrammingQ2(q1n', t', q2n', t', t', t', 0, 7, 1);
            g = interp1(T, G, t');
            testCase.verifyEqual(gam(:), (g(:) - g(1)) / (g(end) - g(1)), 'AbsTol', 1e-14);
        end
    end
end

function d = alignError(testCase, gam)
qw = warp_q_gamma(testCase.q2, gam, testCase.t);
d = sqrt(trapz(testCase.t, (testCase.q1 - qw(:)).^2));
end
