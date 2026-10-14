function [amp, phase] = elastic_depth(f, time, lambda, parallel, nworkers)
% ELASTIC_DISTANCE Calculates the elastic depth between
% functions
% -------------------------------------------------------------------------
% This functions calculates the depth between a set of functions
%
% Usage: [dy, dx] = elastic_depth(f1, f2, time)
%        [dy, dx] = elastic_depth(f1, f2, time, lambda)
%        [dy, dx] = elastic_depth(f1, f2, time, lambda, parallel, nworkers)
%
% Input:
% f: matrix of M time points of N functions, e.g., MxN
% time: vector of size \eqn{N} describing the sample points
% lambda controls amount of warping (default = 0)
% parallel: compute in parallel (default = 1)
% nworkers: number of workers for the parallel pool (0: prompt)
%
% Output
% amp: amplitude depth
% phase: phase depth
arguments
    f double
    time double
    lambda = 0;
    parallel = 1;
    nworkers = 0;
end

M = parfor_workers(parallel, nworkers);

N = size(f,2);

% compute the pairwise distance
phs_dist = zeros(N,N);
amp_dist = zeros(N,N);
parfor (i = 1:N, M)
    for j = i:N
        [da, dp] = elastic_distance(f(:,i),f(:,j),time,lambda);
        
        % y-distance
        amp_dist(i,j) = da;
        
        % x-distance
        phs_dist(i,j) = dp;
    end
end

amp_dist = amp_dist + amp_dist.';
phs_dist = phs_dist + phs_dist.';

amp = 1 ./ (1 + median(amp_dist, 1));
phase = 1 ./ (1 + median(phs_dist, 1));
phase = ((2+pi)/pi) .* (phase - 2/(2+pi));

if parallel == 1
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
end
