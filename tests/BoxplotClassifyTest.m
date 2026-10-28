classdef BoxplotClassifyTest < matlab.unittest.TestCase
    % BOXPLOTCLASSIFYTEST Boxplot region edges and classify of new functions
    %   construct_boxplot stores the region edges it uses, and classify
    %   places new functions in those regions after aligning them to the
    %   fitted template.

    properties
        f
        t
        out
        amp
        ph
    end

    methods (TestClassSetup)
        function setupData(testCase)
            setup_paths;
            S = load('simu_data.mat');
            testCase.f = S.f;
            testCase.t = S.t(:);
            out = fdawarp(testCase.f, testCase.t);
            testCase.out = out.time_warping_median(0, MaxItr=1);
            amp = ampbox(testCase.out);
            testCase.amp = amp.construct_boxplot(.05, 1);
            ph = phbox(testCase.out);
            testCase.ph = ph.construct_boxplot(.05, 1);
        end
    end

    methods (Test)
        function testRegionEdges(testCase)
            for B = {testCase.amp, testCase.ph}
                b = B{1};
                N = numel(b.dist);
                testCase.verifyEqual(b.r50, max(mink(b.dist, ceil(N/2))));
                testCase.verifyEqual(b.ralpha, max(mink(b.dist, round(N*(1-b.alpha)))));
                testCase.verifyTrue(all(b.dist(b.outlier_index) > b.whisker_dis), ...
                    'outliers must lie beyond the whisker distance');
            end
        end

        function testRegionLevels(testCase)
            % the whisker may fall inside ralpha; beyond it is always an outlier
            [level, region] = boxplot_region([0.5 1.5 2.5 4], [1 3 2], .05);
            testCase.verifyEqual(level, [1 2 4 4]);
            testCase.verifyEqual(region(2), "95% central region");
        end

        function testClassifyOutputs(testCase)
            fnew = testCase.f(:, 1:3);
            for B = {testCase.amp, testCase.ph}
                c = B{1}.classify(fnew);
                testCase.verifySize(c.dist, [1 3]);
                testCase.verifyTrue(all(isfinite(c.dist)));
                testCase.verifyTrue(all(ismember(c.level, 1:4)));
                testCase.verifyTrue(all(ismember(c.side, ["Q1" "Q3"])));
                testCase.verifySize(c.gam, size(fnew));
                testCase.verifyLessThan(max(abs(c.gam(1,:))), 1e-12);
                testCase.verifyLessThan(max(abs(c.gam(end,:) - 1)), 1e-12);
            end
        end

        function testAmplitudeMedianAndOutlier(testCase)
            c = testCase.amp.classify(testCase.out.fmean);
            testCase.verifyEqual(c.level, 1, ...
                'the template itself should be in the 50% region');
            c = testCase.amp.classify(10*testCase.f(:, 1));
            testCase.verifyEqual(c.level, 4, ...
                'a 10x scaled function should be an amplitude outlier');
        end
    end
end
