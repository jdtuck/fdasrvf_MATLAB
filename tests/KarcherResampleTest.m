classdef KarcherResampleTest < matlab.unittest.TestCase
    % KARCHERRESAMPLETEST Tests for Karcher means, resampling and distances
    %   Covers SqrtMean, SqrtMeanInverse, ReSampleCurve and elastic_distance.

    properties
        M
        gamId
        f
        t
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            testCase.M = 101;
            testCase.gamId = linspace(0, 1, testCase.M)';
            S = load('simu_data.mat');   % provides f (M x N) and t
            testCase.f = S.f;
            testCase.t = S.t;
        end
    end

    methods (Test)
        function testSqrtMeanOfIdentities(testCase)
            % Karcher mean of identical identity warps is the identity warp
            gam = repmat(testCase.gamId, 1, 5);
            [mu, gam_mu, ~, vec] = SqrtMean(gam);
            testCase.verifyLessThan(max(abs(gam_mu(:) - testCase.gamId(:))), 1e-6, ...
                'SqrtMean of identity warps should be the identity');
            testCase.verifyLessThan(max(abs(mu(:) - 1)), 1e-6, ...
                'Mean psi of identity warps should be all ones');
            testCase.verifyLessThan(max(abs(vec(:))), 1e-6, ...
                'Shooting vectors of identity warps should be ~0');
        end

        function testSqrtMeanInverseOfIdentities(testCase)
            % Inverse Karcher mean of identity warps is the identity warp
            gam = repmat(testCase.gamId, 1, 5);
            gamI = SqrtMeanInverse(gam);
            testCase.verifyLessThan(max(abs(gamI(:) - testCase.gamId(:))), 1e-6, ...
                'SqrtMeanInverse of identity warps should be the identity');
        end

        function testReSampleCurveSize(testCase)
            % Resampling an open curve yields the requested number of points
            T = 200;
            s = linspace(0, 1, T);
            X = [cos(pi*s); sin(pi*s)];
            Xn = ReSampleCurve(X, 50, false);
            testCase.verifyEqual(size(Xn), [2, 50], ...
                'ReSampleCurve produced the wrong size');
            testCase.verifyTrue(all(isfinite(Xn(:))), ...
                'ReSampleCurve produced non-finite values');
        end

        function testReSampleCurveClosedSize(testCase)
            % The closed-curve branch also yields the requested size
            T = 200;
            s = linspace(0, 1, T);
            X = [cos(2*pi*s); sin(2*pi*s)];
            Xn = ReSampleCurve(X, 50, true);
            testCase.verifyEqual(size(Xn), [2, 50], ...
                'ReSampleCurve (closed) produced the wrong size');
            testCase.verifyTrue(all(isfinite(Xn(:))), ...
                'ReSampleCurve (closed) produced non-finite values');
        end

        function testElasticDistanceSelfZero(testCase)
            % A function aligned to itself has ~zero amplitude and phase distance
            [dy, dx] = elastic_distance(testCase.f(:, 1), testCase.f(:, 1), ...
                testCase.t, 0, 'DP1');
            testCase.verifyLessThan(dy, 1e-12, 'Self amplitude distance not ~0');
            testCase.verifyLessThan(dx, 1e-7, 'Self phase distance not ~0');
        end

        function testElasticDistanceNonNegative(testCase)
            % Distances between two distinct functions are finite and non-negative
            [dy, dx] = elastic_distance(testCase.f(:, 1), testCase.f(:, 2), ...
                testCase.t, 0, 'DP1');
            testCase.verifyGreaterThanOrEqual(dy, 0, 'Amplitude distance is negative');
            testCase.verifyGreaterThanOrEqual(dx, 0, 'Phase distance is negative');
            testCase.verifyTrue(isfinite(dy) && isfinite(dx), ...
                'Distances should be finite');
        end

        function testElasticDistanceMethods(testCase)
            % Each optimization method yields finite, non-negative distances
            methods_ = {'DP', 'DP1', 'RBFGS'};
            for k = 1:numel(methods_)
                [dy, dx] = elastic_distance(testCase.f(:, 1), testCase.f(:, 2), ...
                    testCase.t, 0, methods_{k});
                testCase.verifyTrue(isfinite(dy) && isfinite(dx), ...
                    sprintf('Non-finite distance for method %s', methods_{k}));
                testCase.verifyGreaterThanOrEqual(dy, 0, ...
                    sprintf('Negative amplitude distance for method %s', methods_{k}));
                testCase.verifyGreaterThanOrEqual(dx, 0, ...
                    sprintf('Negative phase distance for method %s', methods_{k}));
            end
        end
    end
end
