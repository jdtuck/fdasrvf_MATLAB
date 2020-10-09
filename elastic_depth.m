function [amp, phase] = elastic_depth(f, time, lambda, parallel)
% ELASTIC_DISTANCE Calculates the elastic depth between
% functions
% -------------------------------------------------------------------------
% This functions calculates the depth between a set of functions
%
% Usage: [dy, dx] = elastic_depth(f1, f2, time)
%        [dy, dx] = elastic_depth(f1, f2, time, lambda)
%
% Input:
% f: matrix of M time points of N functions, e.g., MxN
% time: vector of size \eqn{N} describing the sample points
% lambda controls amount of warping (default = 0)
%
% Output
% amp: amplitude depth
% phase: phase depth
if nargin < 3
    lambda = 0;
    parallel = 1;
elseif nargin < 4
    parallel = 1;
end

if parallel == 1
    if isempty(gcp('nocreate'))
        % prompt user for number threads to use
        nThreads = input('Enter number of threads to use: ');
        if nThreads > 1
            parpool(nThreads);
        elseif nThreads > 12 % check if the maximum allowable number of threads is exceeded
            while (nThreads > 12) % wait until user figures it out
                fprintf('Maximum number of threads allowed is 12\n Enter a number between 1 and 12\n');
                nThreads = input('Enter number of threads to use: ');
            end
            if nThreads > 1
                parpool(nThreads);
            end
        end
    end
end

N = size(f,2);

% compute the pairwise distance
phs_dist = zeros(N,N);
amp_dist = zeros(N,N);
parfor i = 1:N
    for j = i:N
        [da, dp] = elastic_distance(f(:,i),f(:,j),time,lambda)
        
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
