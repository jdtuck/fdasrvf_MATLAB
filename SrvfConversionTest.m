classdef SrvfConversionTest < matlab.unittest.TestCase
    % SRVFCONVERSIONTEST Tests for SRVF <-> function and curve <-> SRVF conversions
    %   Covers f_to_srvf, srvf_to_f, curve_to_q and q_to_curve.

    methods (TestClassSetup)
        function setupPaths(~)
            setup_paths;
        end
    end

    methods (Test)
        function testFtoSrvfRoundTrip(testCase)
            % f -> q -> f recovers the original (smooth=false, gradient path)
            t = linspace(0, 1, 101)';
            f = sin(2*pi*t);
            q = f_to_srvf(f, t, false);
            frec = srvf_to_f(q, t, f(1));
            testCase.verifyLessThan(max(abs(frec(:) - f(:))), 1e-6, ...
                'f_to_srvf/srvf_to_f round-trip failed');
        end

        function testFtoSrvfSmoothFalseVsGradient(testCase)
            % smooth=false branch must match the manual SRVF formula exactly
            t = linspace(0, 1, 101)';
            f = sin(2*pi*t) + 0.5*t.^2;
            q = f_to_srvf(f, t, false);
            binsize = mean(diff(t));
            fy = gradient(f, binsize);
            qexp = fy ./ sqrt(abs(fy) + eps);
            testCase.verifyLessThan(max(abs(q(:) - qexp(:))), 1e-12, ...
                'f_to_srvf smooth=false does not match gradient formula');
        end

        function testFtoSrvfConstantIsZero(testCase)
            % A constant function has zero velocity, hence zero SRVF
            t = linspace(0, 1, 101)';
            f = 3.0 * ones(size(t));
            q = f_to_srvf(f, t, false);
            testCase.verifyLessThan(max(abs(q(:))), 1e-6, ...
                'SRVF of a constant function should be ~0');
        end

        function testCurveToQSizeAndScale(testCase)
            % SRVF of a curve is n-by-T and scaled to unit length
            T = 200;
            s = linspace(0, 1, T);
            p = [cos(2*pi*s); sin(2*pi*s)];
            q = curve_to_q(p, false, true);
            testCase.verifyEqual(size(q), [2, T], ...
                'curve_to_q output has wrong size');
            testCase.verifyEqual(sqrt(InnerProd_Q(q, q)), 1, 'AbsTol', 1e-6, ...
                'Scaled SRVF should have unit length');
        end

        function testQtoCurveSizeAndBase(testCase)
            % q_to_curve returns an n-by-T, finite curve starting near origin.
            % NOTE: q_to_curve is not the exact inverse of curve_to_q, so we
            % only check structural properties.
            T = 200;
            s = linspace(0, 1, T);
            p = [cos(2*pi*s); sin(2*pi*s)];
            q = curve_to_q(p, false, true);
            p2 = q_to_curve(q);
            testCase.verifyEqual(size(p2), [2, T], ...
                'q_to_curve output has wrong size');
            testCase.verifyTrue(all(isfinite(p2(:))), ...
                'q_to_curve produced non-finite values');
            testCase.verifyLessThan(max(abs(p2(:, 1))), 1e-9, ...
                'q_to_curve should start at the origin');
        end

        function testCurveToQClosedProjection(testCase)
            % The closed-curve projection branch runs and stays finite
            T = 200;
            s = linspace(0, 1, T);
            p = [cos(2*pi*s); sin(2*pi*s)];
            q = curve_to_q(p, true, true);
            testCase.verifyEqual(size(q), [2, T], ...
                'closed curve_to_q output has wrong size');
            testCase.verifyTrue(all(isfinite(q(:))), ...
                'closed curve_to_q produced non-finite values');
        end
    end
end
