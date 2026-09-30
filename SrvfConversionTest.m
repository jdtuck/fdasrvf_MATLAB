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
            % f -> q -> f recovers the original (smooth=false, spline path)
            t = linspace(0, 1, 101)';
            f = sin(2*pi*t);
            q = f_to_srvf(f, t, false);
            frec = srvf_to_f(q, t, f(1));
            testCase.verifyLessThan(max(abs(frec(:) - f(:))), 1e-6, ...
                'f_to_srvf/srvf_to_f round-trip failed');
        end

        function testFtoSrvfRoundTripComplex(testCase)
            % Oscillatory, sharp and multi-column functions round-trip
            % accurately (errors relative to the range of f)
            t = linspace(0, 1, 401)';
            f = [sin(2*pi*10*t), sin(2*pi*30*t) + 0.3*t, ...
                exp(-200*(t-0.5).^2), t.^3.*sin(40*t)];
            tol = [1e-4, 2e-3, 1e-5, 1e-4];
            q = f_to_srvf(f, t, false);
            frec = srvf_to_f(q, t, f(1,:));
            for k = 1:size(f,2)
                err = max(abs(frec(:,k) - f(:,k))) / (max(f(:,k)) - min(f(:,k)));
                testCase.verifyLessThan(err, tol(k), ...
                    sprintf('round-trip error too large for column %d', k));
            end
        end

        function testFtoSrvfRoundTripConvergence(testCase)
            % Refining the grid reduces the round-trip error (>= 3rd order)
            err = zeros(1,2);
            for j = 1:2
                t = linspace(0, 1, 100*2^(j-1)+1)';
                f = sin(2*pi*5*t);
                q = f_to_srvf(f, t, false);
                err(j) = max(abs(srvf_to_f(q, t, f(1)) - f));
            end
            testCase.verifyGreaterThan(err(1)/err(2), 8, ...
                'round-trip error does not converge fast enough');
        end

        function testFtoSrvfSmoothFalseMatchesSplineDerivative(testCase)
            % smooth=false must match the interpolating spline derivative
            t = linspace(0, 1, 101)';
            f = sin(2*pi*t) + 0.5*t.^2;
            q = f_to_srvf(f, t, false);
            fy = fnval(fnder(csapi(t, f)), t);
            qexp = fy ./ sqrt(abs(fy) + eps);
            testCase.verifyLessThan(max(abs(q(:) - qexp(:))), 1e-12, ...
                'f_to_srvf smooth=false does not match spline derivative');
            % and approximates the true SRSF of f
            fyt = 2*pi*cos(2*pi*t) + t;
            qt = fyt ./ sqrt(abs(fyt) + eps);
            testCase.verifyLessThan(max(abs(q(:) - qt)), 1e-2, ...
                'f_to_srvf smooth=false inaccurate versus analytic SRSF');
        end

        function testFtoSrvfDefaultIsInterpolating(testCase)
            % The default uses the interpolating spline (round-trip safe);
            % smooth=true remains available and differs on noisy data
            t = linspace(0, 1, 101)';
            f = sin(2*pi*t) + 0.05*sin(2*pi*40*t);
            testCase.verifyEqual(f_to_srvf(f, t), f_to_srvf(f, t, false));
            testCase.verifyGreaterThan(max(abs(f_to_srvf(f, t) - f_to_srvf(f, t, true))), 1e-3);
        end

        function testSrvfToFStartsAtFo(testCase)
            % srvf_to_f honors the initial value of each function
            t = linspace(0, 1, 101)';
            q = [ones(101,1), -ones(101,1)];
            f = srvf_to_f(q, t, [2, -3]);
            testCase.verifyEqual(f(1,:), [2, -3], 'AbsTol', 1e-12);
            testCase.verifyEqual(f(end,:), [3, -4], 'AbsTol', 1e-10);
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
