function [IndicatorMatrix, curvatures, heights, heights2] = PreprocessingForPPD(t, lam, Labels, Locs, labelMax, FNm, th)

K = length(lam);

% Initialize output matrices
IndicatorMatrix = nan(K, labelMax);
curvatures      = zeros(K, labelMax);
heights         = nan(K, labelMax);

% Assume t is uniformly spaced; compute time step
dx = t(2) - t(1);

% Process each function
for i = 1:K
    % Extract the function values for the current parameter lambda
    fnm = FNm(:, i);

    % Compute negative curvature (second derivative)
    negCurvature = -gradient(gradient(fnm, dx), dx);

    % Ensure non-negative curvature values
    negCurvature = max(0, negCurvature);

    % Normalize negative curvature to [0, 1] if possible
    maxNegCurvature = max(negCurvature);
    if maxNegCurvature > 0
        negCurvature = negCurvature / maxNegCurvature;
    end

    % Retrieve peak locations and labels for the current function
    locsCurrent   = Locs{i};
    labelsCurrent = Labels{i};

    % Select negative curvature values at specified peak locations
    negCurvSelected = negCurvature(locsCurrent);

    % Update curvatures and heights matrices at the appropriate labels
    curvatures(i, labelsCurrent) = negCurvSelected;
    heights(i, labelsCurrent)    = fnm(locsCurrent);

    % Apply threshold to select significant peaks based on curvature
    significantLabels = labelsCurrent(negCurvSelected >= th);

    % Update the indicator matrix for significant peaks
    IndicatorMatrix(i, significantLabels) = 1;
end

% Compute heights2 by multiplying heights with the indicator matrix
heights2 = IndicatorMatrix .* heights;

% Replace zeros in the indicator matrix with NaN for clarity in plotting
IndicatorMatrix(isnan(IndicatorMatrix)) = NaN;

end