function M = parfor_workers(parallel, nworkers)
% PARFOR_WORKERS Number of workers to pass to a parfor loop
% -------------------------------------------------------------------------
% Usage:  M = parfor_workers(parallel)
%         M = parfor_workers(parallel, nworkers)
%
% Use as parfor (i = 1:n, M) so the loop runs in parallel only when
% requested; M = 0 runs the loop serially on the client.
%
% Input:
% parallel: run in parallel (starts a pool if one is not already running)
% nworkers: number of workers for the parallel pool (0: prompt)
%
% Output:
% M: maximum number of workers for the parfor loop

if nargin < 2
    nworkers = 0;
end

M = 0;
if parallel
    start_parpool(nworkers);
    pool = gcp('nocreate');
    if ~isempty(pool)
        M = pool.NumWorkers;
    end
end
end
