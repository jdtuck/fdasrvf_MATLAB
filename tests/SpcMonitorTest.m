classdef SpcMonitorTest < matlab.unittest.TestCase
    % SPCMONITORTEST Tests for elastic functional SPC (fdaspc)
    %   Checks the T^2/SPE statistics and limits for amplitude and phase,
    %   and that amplitude and phase faults raise alarms.

    properties
        warp
        t
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            d = load('simu_data.mat');
            testCase.t = d.t(:);
            obj = fdawarp(d.f, testCase.t);
            testCase.warp = obj.time_warping(0);
        end
    end

    methods (Test)
        function testInSampleT2MatchesCoef(testCase)
            % without re-alignment T^2 is the scaled fPCA coefficient norm
            s = fdaspc(testCase.warp).train(realign=false);
            testCase.verifyEqual(s.amp.T2, ...
                sum(s.vpca.coef.^2 ./ s.vpca.latent.', 2), 'AbsTol', 1e-10);
            testCase.verifyEqual(s.ph.T2, ...
                sum(s.hpca.coef.^2 ./ s.hpca.latent.', 2), 'AbsTol', 1e-10);
        end

        function testContributionsSum(testCase)
            s = fdaspc(testCase.warp).train(no=[3 2]);
            [~, out] = s.monitor(testCase.warp.f(:, 1:4));
            for m = {out.amp, out.ph}
                testCase.verifyEqual(sum(m{1}.contT2, 2), m{1}.T2, 'RelTol', 1e-8);
                testCase.verifyEqual(sum(m{1}.contSPE, 2), m{1}.SPE, 'RelTol', 1e-8);
            end
        end

        function testTrainingMostlyInControl(testCase)
            s = fdaspc(testCase.warp).train(alpha=0.05);
            [~, out] = s.monitor(testCase.warp.f);
            testCase.verifyLessThan(mean(out.alarm), 0.5);
        end

        function testFaultsAlarm(testCase)
            s = fdaspc(testCase.warp).train(alpha=0.05);
            f = testCase.warp.f(:, 1:3);
            t = testCase.t;
            g = min(max(t + 0.08*sin(2*pi*t), 0), 1);
            fp = zeros(size(f));
            for i = 1:3
                fp(:, i) = warp_f_gamma(f(:, i), g, t);
            end
            [s, out] = s.monitor([1.5*f fp]);
            testCase.verifyTrue(all(out.amp.alarm(1:3)), 'amplitude fault missed');
            testCase.verifyTrue(all(out.ph.alarm(4:6)), 'phase fault missed');
            testCase.verifyEqual(s.mon.alarm, out.alarm);
        end
    end
end
