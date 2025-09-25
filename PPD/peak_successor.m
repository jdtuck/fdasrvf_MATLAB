function [labels2, labelMax] = peak_successor(f1, f2, labels1, labelMax, smooth_parameter)

if isempty(f1) || isempty(f2)
    error('Input data f1 and f2 must be non-empty.');
end

f1 = smoothdata(f1,'movmean',smooth_parameter);
f2 = smoothdata(f2,'movmean',smooth_parameter);

% Combine f1 and f2 into a 3D array and compute the mean across the second dimension
F = cat(3, f1, f2);
fm = squeeze(mean(F, 2));  % Resulting in an n x 2 matrix

fm = fillmissing(fm, 'linear');

% Compute peak ranges and labels for fm(:,1)
[ranges, idx_max1] = computePeakRanges(fm(:, 1));

% Compute indices of local maxima in fm(:,2)
idx_max2 = find(islocalmax(fm(:, 2)));

if isempty(idx_max1)
    labels2 = labelMax + (1:numel(idx_max2));

else
    % Assign labels to peaks in fm(:,2) based on matching ranges in fm(:,1)
    labels2 = assignLabelsToPeaks(idx_max2, ranges, idx_max1, labels1);

    % Ensure no overlapping labels in labels2
    labels2 = resolveOverlappingLabels(labels2, idx_max1, idx_max2, labels1);

    % Assign new labels to unmatched peaks in fm(:,2)
    unmatched = labels2 == 0;
    labels2(unmatched) = labelMax + (1:sum(unmatched));
    labelMax = labelMax + sum(unmatched);

end

% --- Nested Helper Functions ---

    function [ranges, idx_max] = computePeakRanges(data)
        % Computes peak ranges defined by adjacent minima in the data
        idx_max = find(islocalmax(data), 1);

        if isempty(idx_max)
            warning('No peaks found in f1.');
            idx_max = [];
            ranges = [];
            return
        end

        idx_max = find(islocalmax(data, 'FlatSelection', 'first'));
        idx_min = unique([1; length(data); find(islocalmin(data))]);
        ranges = arrayfun(@(idx) [max(idx_min(idx_min < idx)), min(idx_min(idx_min > idx))], ...
            idx_max, 'UniformOutput', false);
        ranges = vertcat(ranges{:});
        ranges(ranges(:, 1) == ranges(:, 2), :) = [];  % Remove degenerate ranges
    end

    function labels = assignLabelsToPeaks(idx_max2, ranges, idx_max1, labels1)
        % Assigns labels to peaks in idx_max2 based on matching ranges in idx_max1
        labels = zeros(1, numel(idx_max2));
        for i = 1:numel(idx_max2)
            % Find the range in fm(:,1) that contains the current peak in fm(:,2)
            in_range = idx_max2(i) >= ranges(:, 1) & idx_max2(i) <= ranges(:, 2);
            matching_ranges = find(in_range);

            if ~isempty(matching_ranges)
                % Choose the closest peak if multiple ranges match
                if numel(matching_ranges) > 1
                    [~, closest_idx] = min(abs(idx_max1(matching_ranges) - idx_max2(i)));
                    matching_range = matching_ranges(closest_idx);
                else
                    matching_range = matching_ranges;
                end
                labels(i) = labels1(matching_range);
            end
        end
    end

    function labels = resolveOverlappingLabels(labels, idx_max1, idx_max2, labels1)
        % Ensures no overlapping labels in the assigned labels
        unique_labels = unique(labels(labels > 0));
        for label = unique_labels
            duplicates = find(labels == label);
            if numel(duplicates) > 1
                % Keep the closest peak and reset others
                distances = abs(idx_max1(labels1 == label) - idx_max2(duplicates));
                [~, min_idx] = min(distances);
                duplicates(min_idx) = [];
                labels(duplicates) = 0;
            end
        end
    end

end