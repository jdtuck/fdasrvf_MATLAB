classdef DynamicProgrammingTest < matlab.unittest.TestCase
    % DYNAMICPROGRAMMINGTEST Accuracy of the DP1 alignment on a known warp
    %   f2 = f1 o gam_true, so aligning q2 to q1 should recover the inverse
    %   of gam_true.  The DP cost interpolates q2 between grid points;
    %   rounding to the nearest grid point instead left the warp up to 0.63
    %   grid cells off and the fit worse than the true warp's.

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
    end
end

function d = alignError(testCase, gam)
qw = warp_q_gamma(testCase.q2, gam, testCase.t);
d = sqrt(trapz(testCase.t, (testCase.q1 - qw(:)).^2));
end
