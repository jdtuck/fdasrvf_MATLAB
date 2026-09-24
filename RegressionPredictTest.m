classdef RegressionPredictTest < matlab.unittest.TestCase
    % REGRESSIONPREDICTTEST Tests for elastic PCR model fitting and prediction
    %   Covers elastic_pcr_regression, elastic_lpcr_regression and
    %   elastic_mlpcr_regression, including prediction on new data, for each
    %   fPCA method.

    properties
        f
        t
    end

    properties (TestParameter)
        pcaMethod = {'vert', 'horiz', 'combined'};
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            setup_paths;
            d = load('data/simu_data.mat');
            testCase.f = d.f;
            testCase.t = d.t;
        end
    end

    methods (Test)
        function testPcrPredict(testCase, pcaMethod)
            N = size(testCase.f, 2);
            y = linspace(-1, 1, N).';
            model = elastic_pcr_regression(testCase.f, y, testCase.t);
            model = model.calc_model(pcaMethod, 0.99, parallel=0);

            model = model.predict();
            yfit = model.y_pred;
            testCase.verifySize(yfit, [N 1]);
            testCase.verifyTrue(all(isfinite(yfit)));

            newdata = new_data(testCase.f(:,1:5), testCase.t, y(1:5));
            model = model.predict(newdata);
            testCase.verifySize(model.y_pred, [5 1]);
            testCase.verifyTrue(all(isfinite(model.y_pred)));
            testCase.verifyTrue(isfinite(model.SSE));
        end

        function testLpcrPredict(testCase, pcaMethod)
            N = size(testCase.f, 2);
            y = ones(N, 1);
            y(1:2:end) = -1;
            model = elastic_lpcr_regression(testCase.f, y, testCase.t);
            model = model.calc_model(pcaMethod, 0.99, parallel=0);

            model = model.predict();
            testCase.verifySize(model.y_labels, [N 1]);
            testCase.verifyTrue(all(ismember(model.y_labels, [-1 1])));
            testCase.verifyGreaterThanOrEqual(model.PC, 0);
            testCase.verifyLessThanOrEqual(model.PC, 1);

            newdata = new_data(testCase.f(:,1:5), testCase.t, y(1:5));
            model = model.predict(newdata);
            testCase.verifySize(model.y_labels, [5 1]);
            testCase.verifyTrue(all(ismember(model.y_labels, [-1 1])));
            testCase.verifyTrue(isfinite(model.PC));
        end

        function testMlpcrPredict(testCase, pcaMethod)
            N = size(testCase.f, 2);
            y = mod((0:N-1).', 3) + 1;
            model = elastic_mlpcr_regression(testCase.f, y, testCase.t);
            model = model.calc_model(pcaMethod, 0.99, parallel=0);
            testCase.verifySize(model.alpha, [1 3]);
            testCase.verifySize(model.b, [size(model.pca.coef,2) 3]);

            model = model.predict();
            testCase.verifySize(model.y_labels, [N 1]);
            testCase.verifyTrue(all(ismember(model.y_labels, 1:3)));
            testCase.verifyGreaterThanOrEqual(model.PC, 0);
            testCase.verifyLessThanOrEqual(model.PC, 1);

            newdata = new_data(testCase.f(:,1:5), testCase.t, y(1:5));
            model = model.predict(newdata);
            testCase.verifySize(model.y_labels, [5 1]);
            testCase.verifyTrue(all(ismember(model.y_labels, 1:3)));
            testCase.verifyTrue(isfinite(model.PC));
        end

        function testMlpcrLabelsFollowLargestScore(testCase)
            % Training fits P(class j) proportional to exp(score j), so the
            % predicted label must be the class with the largest score
            N = size(testCase.f, 2);
            y = mod((0:N-1).', 3) + 1;
            model = elastic_mlpcr_regression(testCase.f, y, testCase.t);
            model = model.calc_model('vert', 0.99, parallel=0);
            model = model.predict();

            scores = model.alpha + model.pca.coef*model.b;
            [~, expected] = max(scores, [], 2);
            testCase.verifyEqual(model.y_labels, expected);
        end
    end
end

function newdata = new_data(f, t, y)
newdata.f = f;
newdata.time = t;
newdata.y = y;
newdata.smooth = 0;
newdata.sparam = 25;
end
