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
            % "DP" integrates the linearly interpolated SRVFs exactly over a
            % 10-point neighborhood; the old step-function cost over 7 points
            % fit 1.5x worse than the true warp's
            gam = optimum_reparam(testCase.q1, testCase.q2, testCase.t, 0, 'DP');
            dy = alignError(testCase, gam(:));
            dyTrue = alignError(testCase, testCase.gamInv);
            testCase.verifyLessThan(dy, 1.1 * dyTrue, ...
                'DP alignment should fit about as well as the true warp');
        end

        function testRBFGSFitsAsWellAsTrueWarp(testCase)
            % "RBFGS" refines the "DP" warp; started from the identity it
            % stopped half a grid cell off
            gam = optimum_reparam(testCase.q1, testCase.q2, testCase.t, 0, 'RBFGS');
            M = length(testCase.t);
            testCase.verifyLessThan(max(abs(gam(:) - testCase.gamInv)) * (M-1), 0.5, ...
                'RBFGS warp should be within half a grid cell of the true warp');
            dy = alignError(testCase, gam(:));
            dyTrue = alignError(testCase, testCase.gamInv);
            testCase.verifyLessThan(dy, 1.1 * dyTrue, ...
                'RBFGS alignment should fit about as well as the true warp');
        end

        function testRBFGSCostMatchesDP(testCase)
            % On random pairs the identity start stopped in local minima, up
            % to 3x the "DP" distance.  Scored on DP's own objective, RBFGS
            % is below DP on most pairs; from the coarse DP start it can land
            % in another minimum, and the worst of these pairs is 1.30x
            t = testCase.t;
            rng(0);
            r = zeros(10, 1);
            for i = 1:10
                f1 = bumps(t);
                f2 = bumps(t);
                q1 = f_to_srvf(f1, t);
                q2 = f_to_srvf(f2, t);
                r(i) = dpCost(q1, q2, optimum_reparam(q1, q2, t, 0, 'RBFGS'), t) / ...
                    dpCost(q1, q2, optimum_reparam(q1, q2, t, 0, 'DP'), t);
            end
            testCase.verifyLessThan(median(r), 1.0, ...
                'RBFGS should fit at least as well as DP on average');
            testCase.verifyLessThan(max(r), 1.35, ...
                'RBFGS should fit within 35% of DP on every pair');
        end

        function testRBFGSSolverOutputs(testCase)
            % c_rlbfgs starts from gam0 when given one and reports its stop
            t = testCase.t;
            q1n = testCase.q1 / norm(testCase.q1);
            q2n = testCase.q2 / norm(testCase.q2);
            [g, info] = c_rlbfgs(q1n, q2n, t, 30, 0, 0, testCase.gamInv);
            testCase.verifyEqual(size(info), [1 4]);
            testCase.verifyTrue(ismember(info(2), [1 2 3 4]), 'unknown stop reason');
            testCase.verifyTrue(all(diff(g) >= 0), 'warp should be monotone');
            testCase.verifyLessThan(max(abs(g(:) - testCase.gamInv)) * (length(t)-1), 0.5, ...
                'started at the true warp, RBFGS should stay near it');
        end

        function testRBFGSStartsAtGam0(testCase)
            % Started from a DP warp, RBFGS starts there and then descends;
            % the old solver started up to 2.9x above its cost and stalled
            [q1, q2, t] = growthPair();
            q1 = q1 / norm(q1);
            q2 = q2 / norm(q2);
            gd = optimum_reparam(q1, q2, t, 0, 'DP');
            t1 = linspace(0, 1, length(t))';
            [g1, info1] = c_rlbfgs(q1, q2, t1, 1, 0, 0, gd(:));
            [~, info30] = c_rlbfgs(q1, q2, t1, 30, 0, 0, gd(:));
            psi1 = gam_to_psi(g1(:), false);
            psid = gam_to_psi(gd(:), false);
            dx = acos(min(trapz(t1, psi1 .* psid), 1));
            testCase.verifyLessThan(dx, 0.15, ...
                'one iteration from the DP warp should stay near it');
            testCase.verifyLessThan(info30(3), info1(3), ...
                'further iterations should lower the cost');
        end

        function testRBFGSMatchesDPOnData(testCase)
            % On the growth velocities, whose sharp spurts made the last
            % solver 1.04x DP's cost in the median, RBFGS is below DP's cost
            % in the median and at most 1.19x on any consecutive pair
            S = load(fullfile(fileparts(which('optimum_reparam')), 'data', ...
                'growth_female_vel.mat'));
            t = S.t(:);
            q = f_to_srvf(S.f, t);
            N = size(q, 2);
            r = zeros(N-1, 1);
            for i = 1:N-1
                gd = optimum_reparam(q(:,i), q(:,i+1), t, 0, 'DP');
                gr = optimum_reparam(q(:,i), q(:,i+1), t, 0, 'RBFGS');
                r(i) = dpCost(q(:,i), q(:,i+1), gr, t) / dpCost(q(:,i), q(:,i+1), gd, t);
            end
            testCase.verifyLessThan(median(r), 1.0, ...
                'RBFGS should fit at least as well as DP on average');
            testCase.verifyLessThan(max(r), 1.25, ...
                'RBFGS should fit within 25% of DP on every pair');
        end

        function testSolverArguments(testCase)
            % scl reaches DynamicProgrammingQ, and "DP" defaults to nbhd_dim = 10
            t = testCase.t;
            q1n = testCase.q1 / norm(testCase.q1);
            q2n = testCase.q2 / norm(testCase.q2);
            gam = optimum_reparam(testCase.q1, testCase.q2, t, 0, 'DP1', 0, 0, 7, 'roughness', 1);
            g = DynamicProgrammingQ(q2n', q1n', 0, 1, 1);
            testCase.verifyEqual(gam(:), (g(:) - g(1)) / (g(end) - g(1)), 'AbsTol', 1e-14);
            gam = optimum_reparam(testCase.q1, testCase.q2, t, 0, 'DP');
            [G, T] = DynamicProgrammingQ2(q1n', t', q2n', t', t', t', 0, 10, 1);
            g = interp1(T, G, t');
            testCase.verifyEqual(gam(:), (g(:) - g(1)) / (g(end) - g(1)), 'AbsTol', 1e-14);
        end
    end
end

function d = alignError(testCase, gam)
qw = warp_q_gamma(testCase.q2, gam, testCase.t);
d = sqrt(trapz(testCase.t, (testCase.q1 - qw(:)).^2));
end

function E = dpCost(q1, q2, gam, t)
% DP's objective at gam: the unit norm, linearly interpolated SRVFs, which
% DynamicProgrammingQ2 integrates exactly, here on a 20x finer grid
M = length(t);
t1 = linspace(0, 1, M)';
tf = linspace(0, 1, (M-1)*20 + 1)';
a = interp1(t1, q1(:), tf); a = a / sqrt(trapz(tf, a.^2));
b = interp1(t1, q2(:), tf); b = b / sqrt(trapz(tf, b.^2));
E = trapz(tf, (a - warp_q_gamma(b, interp1(t1, gam(:), tf), tf)').^2);
end

function [q1, q2, t] = growthPair()
% a growth velocity pair on which the old gam0 start stalled
S = load(fullfile(fileparts(which('optimum_reparam')), 'data', 'growth_female_vel.mat'));
t = S.t(:);
q = f_to_srvf(S.f, t);
q1 = q(:, 11);
q2 = q(:, 12);
end

function f = bumps(t)
% three random Gaussian bumps
f = zeros(size(t));
for j = 1:3
    f = f + randn * exp(-(t - rand).^2 / (2*(0.05 + 0.1*rand)^2));
end
end
