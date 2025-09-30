function Clt2 = getPersistentPeaks(IndicatorMatrix)

% Validate input
if isempty(IndicatorMatrix)
    Clt2 = [];
    return;
end

% Count the number of ones (occurrences) for each peak (ignore NaNs)
occurrenceCounts = sum(IndicatorMatrix, 2, "omitnan");


data = [occurrenceCounts; 0];
if all(data == 0) || isscalar(unique(occurrenceCounts))
    Clt2 = [];
    return;
end

% Compute pairwise distances between observations
pairwiseDistances = pdist(data, 'euclidean');
Y = linkage(pairwiseDistances, 'ward');

% Cluster the data into the specified number of clusters
clusterAssignments = cluster(Y, 'maxclust', 2);
referenceCluster = clusterAssignments(end);
clusterAssignments(end) = [];

% Identify indices where cluster assignments differ from the reference
Clt2 = find(clusterAssignments ~= referenceCluster)';

end