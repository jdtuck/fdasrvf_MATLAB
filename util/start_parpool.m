function start_parpool(nworkers)
% START_PARPOOL Start a parallel pool if one is not already running
% -------------------------------------------------------------------------
% Usage:  start_parpool()
%         start_parpool(nworkers)
%
% Input:
% nworkers: number of workers to use. If empty or 0 (default), the user is
%           prompted for the number of workers. Values larger than the
%           number of available cores are capped.

if nargin < 1
    nworkers = 0;
end

if ~isempty(gcp('nocreate'))
    return
end

if isempty(nworkers) || nworkers < 1
    nworkers = input('Enter number of threads to use: ');
end

maxWorkers = feature('numcores');
if nworkers > maxWorkers
    fprintf('Requested %d workers, but only %d cores available. Using %d.\n', ...
        nworkers, maxWorkers, maxWorkers);
    nworkers = maxWorkers;
end

if nworkers > 1
    parpool(nworkers);
end
end
