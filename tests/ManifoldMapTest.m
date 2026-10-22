classdef ManifoldMapTest < matlab.unittest.TestCase
    % MANIFOLDMAPTEST Tests for Hilbert-sphere maps and functional products
    %   Covers exp_map, inv_exp_map, inner_product and L2norm.

    properties
        M
        psi
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            setup_paths;
            testCase.M = 101;
            testCase.psi = ones(testCase.M, 1);   % identity-warp SRVF point
        end
    end

    methods (Test)
        function testInvExpMapSelfZero(testCase)
            % The inverse-exp map of a point with itself is the zero tangent
            out = inv_exp_map(testCase.psi, testCase.psi);
            testCase.verifyLessThan(max(abs(out(:))), 1e-12, ...
                'inv_exp_map(psi,psi) should be zero');
        end

        function testExpMapZeroReturnsBase(testCase)
            % The exp map of a zero tangent returns the base point
            out = exp_map(testCase.psi, zeros(testCase.M, 1));
            testCase.verifyLessThan(max(abs(out(:) - testCase.psi(:))), 1e-12, ...
                'exp_map(psi,0) should return psi');
        end

        function testExpInvExpInverse(testCase)
            % inv_exp_map o exp_map recovers a small tangent vector.
            % v must be tangent (orthogonal) to psi under inner_product.
            v0 = sin(2*pi*linspace(0, 1, testCase.M))';
            ip = inner_product(v0, testCase.psi);
            v = v0 - ip * testCase.psi;      % project onto tangent space
            v = 0.1 * v / L2norm(v);         % small step
            p2 = exp_map(testCase.psi, v);
            vrec = inv_exp_map(testCase.psi, p2);
            testCase.verifyLessThan(max(abs(vrec(:) - v(:))), 1e-6, ...
                'inv_exp_map o exp_map did not recover the tangent vector');
        end

        function testL2normOnes(testCase)
            % L2 norm of the constant one function over [0,1] is 1
            testCase.verifyEqual(L2norm(ones(testCase.M, 1)), 1, 'AbsTol', 1e-12, ...
                'L2norm(ones) should be 1');
        end

        function testL2normZeros(testCase)
            testCase.verifyEqual(L2norm(zeros(testCase.M, 1)), 0, 'AbsTol', 1e-15, ...
                'L2norm(zeros) should be 0');
        end

        function testInnerProductOnes(testCase)
            testCase.verifyEqual(inner_product(ones(testCase.M, 1), ...
                ones(testCase.M, 1)), 1, 'AbsTol', 1e-12, ...
                'inner_product(ones,ones) should be 1');
        end

        function testInnerProductSymmetry(testCase)
            a = sin(2*pi*linspace(0, 1, testCase.M))';
            b = cos(2*pi*linspace(0, 1, testCase.M))';
            testCase.verifyEqual(inner_product(a, b), inner_product(b, a), ...
                'AbsTol', 1e-14, 'inner_product should be symmetric');
        end
    end
end
