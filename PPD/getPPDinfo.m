function [IndicatorMatrix, Curvatures, Heights, Locs, Labels, FNm] = getPPDinfo(t, Fa, lam, th)

n_lams = length(lam);

% Compute the mean of each function in FN_temp over its rows
FNm = cell2mat(cellfun(@(x) mean(x,2), Fa, 'UniformOutput', false)');

% Find indices of local maxima in the first function's mean
idxMaxFirst = find(islocalmax(FNm(:, 1)));

% Initialize Labels and Locations for the first function
Labels = cell(1, n_lams);
Locs = cell(1, n_lams);
Labels{1} = 1:length(idxMaxFirst);
Locs{1} = idxMaxFirst;

% Initialize the maximum label number
labelMax = max(Labels{1});

% Process each function to assign labels and locate peaks
for i = 1:n_lams - 1

    currentLabel = Labels{i};
    [Labels{i + 1}, labelMax] = peak_successor(Fa{i}, Fa{i + 1}, currentLabel, labelMax, 1);

    % Find peak locations in the next function's mean
    FNmNextMean = mean(Fa{i + 1}, 2);
    idxMaxNext = find(islocalmax(FNmNextMean));
    Locs{i + 1} = idxMaxNext';

    % Update the mean function matrix
    FNm(:, i+1) = FNmNextMean;
end

% Preprocess data to compute IndicatorMatrix, Curvatures, and Heights
[IndicatorMatrix, Curvatures, Heights, Heights2] = PreprocessingForPPD(t, lam, Labels, Locs, labelMax, FNm, th);

% Check if all peaks are ignored based on the threshold
if all(isnan(Heights2), 'all')
    warning('All peaks are ignored. A smaller threshold is required.');
    return;
end

end