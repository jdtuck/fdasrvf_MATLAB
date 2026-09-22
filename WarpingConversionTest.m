classdef WarpingConversionTest < matlab.unittest.TestCase
    % WARPINGCONVERSIONTEST Tests for warping-coordinate conversions
    %   Covers gam_to_psi/psi_to_gam, gam_to_v/v_to_gam and gam_to_h/h_to_gam.
    %   All use smooth=false so the numerics are exact and no Curve Fitting
    %   Toolbox 'fit' call is exercised.

    properties
        M
        gamId
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            setup_paths;
            testCase.M = 101;
            testCase.gamId = linspace(0, 1, testCase.M)';
        end
    end

    methods (Test)
        function testPsiIdentityIsOnes(testCase)
            % The identity warp maps to psi == 1 everywhere
            psi = gam_to_psi(testCase.gamId, false);
            testCase.verifyLessThan(max(abs(psi(:) - 1)), 1e-6, ...
                'psi of the identity warp should be all ones');
        end

        function testGamToPsiToGam(testCase)
            % gam -> psi -> gam recovers the identity warp
            psi = gam_to_psi(testCase.gamId, false);
            gam = psi_to_gam(psi);
            testCase.verifyLessThan(max(abs(gam(:) - testCase.gamId(:))), 1e-6, ...
                'gam_to_psi/psi_to_gam round-trip failed');
        end

        function testGamToVToGam(testCase)
            % Identity warp has zero shooting vector and round-trips
            v = gam_to_v(testCase.gamId, false, false);
            testCase.verifyLessThan(max(abs(v(:))), 1e-6, ...
                'Shooting vector of the identity warp should be ~0');
            gam = v_to_gam(v);
            testCase.verifyLessThan(max(abs(gam(:) - testCase.gamId(:))), 1e-6, ...
                'gam_to_v/v_to_gam round-trip failed');
        end

        function testGamToHToGam(testCase)
            % Identity warp has zero h coordinate and round-trips
            h = gam_to_h(testCase.gamId, false);
            testCase.verifyLessThan(max(abs(h(:))), 1e-6, ...
                'h coordinate of the identity warp should be ~0');
            gam = h_to_gam(h);
            testCase.verifyLessThan(max(abs(gam(:) - testCase.gamId(:))), 1e-6, ...
                'gam_to_h/h_to_gam round-trip failed');
        end

        function testNonIdentityGamRoundTripPsi(testCase)
            % A non-trivial monotone warp round-trips through psi, confirming
            % the round-trip tests are not merely identity-preserving.
            s = linspace(0, 1, testCase.M)';
            gam0 = (exp(3*s) - 1) / (exp(3) - 1);
            psi = gam_to_psi(gam0, false);
            gam = psi_to_gam(psi);
            testCase.verifyLessThan(max(abs(gam(:) - gam0(:))), 1e-4, ...
                'Non-identity warp did not round-trip through psi');
            testCase.verifyGreaterThan(max(abs(gam0(:) - s)), 1e-2, ...
                'Test warp is too close to identity to be meaningful');
        end
    end
end
