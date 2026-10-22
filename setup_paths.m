% setup paths
root = fileparts(mfilename('fullpath'));
subdirs = {'PNS','PPD','basis','data','manifold','mex','util'};
dirs = {};
for k = 1:numel(subdirs)
    dirs = [dirs, strsplit(genpath(fullfile(root,subdirs{k})), pathsep)]; %#ok<AGROW>
end
dirs = [dirs, {root}];

% only add folders not already on the path; re-adding folders that the
% fdasrvf package/project already put on the path reorders them and
% triggers MATLAB's package folder ordering warning
dirs = dirs(~cellfun('isempty',dirs) & ~ismember(dirs, strsplit(path, pathsep)));
if ~isempty(dirs)
    addpath(dirs{:})
end
clear root subdirs dirs k
