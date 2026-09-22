classdef WarpApplyTest < matlab.unittest.TestCase
    % WARPAPPLYTEST Tests for applying warps and for gamma inversion
    %   Covers warp_f_gamma, warp_q_gamma, warp_srvf_gamma and invertGamma.

    properties
        M
        t
        gamId
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            setup_paths;
            testCase.M = 101;
            testCase.t = linspace(0, 1, testCase.M);
            testCase.gamId = linspace(0, 1, testCase.M);
        end
    end

    methods (Test)
        function testWarpFGammaIdentity(testCase)
            % Warping a function by the identity gamma is a no-op
            f = sin(2*pi*testCase.t);
            fw = warp_f_gamma(f, testCase.gamId, testCase.t);
            testCase.verifyLessThan(max(abs(fw(:) - f(:))), 1e-6, ...
                'warp_f_gamma with identity gamma is not a no-op');
        end

        function testWarpQGammaIdentity(testCase)
            % Warping an SRVF by the identity gamma is a no-op (gam_dev == 1)
            q = sin(2*pi*testCase.t);
            qw = warp_q_gamma(q, testCase.gamId, testCase.t);
            testCase.verifyLessThan(max(abs(qw(:) - q(:))), 1e-6, ...
                'warp_q_gamma with identity gamma is not a no-op');
        end

        function testWarpSrvfGammaIdentityNoScale(testCase)
            % scale=false leaves the SRVF unchanged under the identity warp
            T = testCase.M;
            s = linspace(0, 1, T);
            q = [cos(2*pi*s); sin(2*pi*s)];
            qw = warp_srvf_gamma(q, testCase.gamId, false);
            testCase.verifyLessThan(max(abs(qw(:) - q(:))), 1e-6, ...
                'warp_srvf_gamma (scale=false) with identity gamma is not a no-op');
        end

        function testWarpSrvfGammaIdentityScaled(testCase)
            % scale=true renormalizes the warped SRVF to unit length
            T = testCase.M;
            s = linspace(0, 1, T);
            q = [cos(2*pi*s); sin(2*pi*s)];
            qw = warp_srvf_gamma(q, testCase.gamId, true);
            testCase.verifyEqual(sqrt(InnerProd_Q(qw, qw)), 1, 'AbsTol', 1e-6, ...
                'warp_srvf_gamma (scale=true) should return unit length');
        end

        function testInvertGammaLinearIsIdentity(testCase)
            % Inverting the identity warp gives the identity warp
            gami = invertGamma(testCase.gamId);
            testCase.verifyLessThan(max(abs(gami(:) - testCase.gamId(:))), 1e-12, ...
                'invertGamma of the identity should be the identity');
        end

        function testInvertGammaDouble(testCase)
            % Inverting twice recovers the original warp
            g = linspace(0, 1, testCase.M).^2;
            gii = invertGamma(invertGamma(g));
            testCase.verifyLessThan(max(abs(gii(:) - g(:))), 1e-4, ...
                'Double gamma inversion did not recover the original');
        end

        function testWarpComposeInverse(testCase)
            % Warping by gamma then by its inverse is ~identity on the interior
            % (endpoints are excluded because interpolation can be unstable there)
            f = sin(2*pi*testCase.t);
            g = linspace(0, 1, testCase.M).^1.5;
            gi = invertGamma(g);
            fw = warp_f_gamma(f, g, testCase.t);
            fwi = warp_f_gamma(fw(:)', gi, testCase.t);
            idx = 5:(testCase.M - 4);
            testCase.verifyLessThan(max(abs(fwi(idx) - f(idx)')), 1e-3, ...
                'warp followed by inverse warp did not recover the function');
        end
    end
end
