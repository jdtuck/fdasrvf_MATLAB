function [level, region] = boxplot_region(dist, edges, alpha)
% BOXPLOT_REGION Boxplot region of a distance to the median
% -------------------------------------------------------------------------
% Usage: [level, region] = boxplot_region(dist, edges, alpha)
%
% Input:
% dist: vector of distances to the median
% edges: [r50, ralpha, whisker_dis] from ampbox/phbox
% alpha: quantile level used for the (1-alpha) region
%
% Output:
% level: 1 = 50% region, 2 = (1-alpha) region, 3 = within whiskers,
%        4 = outlier
% region: level as text
names = ["50% central region", ...
    sprintf("%g%% central region", 100*(1-alpha)), ...
    "within whiskers", "outlier"];
% the whisker can fall inside ralpha, so test it first: construct_boxplot
% flags every distance beyond it as an outlier
dist = dist(:)';
level = 3*ones(size(dist));
level(dist <= edges(2)) = 2;
level(dist <= edges(1)) = 1;
level(dist > edges(3)) = 4;
region = names(level);
end
