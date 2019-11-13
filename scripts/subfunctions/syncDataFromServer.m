function dataDirGitRepo = syncDataFromServer(projectPath)
% Function to rsync data from server to git repo, since we don't want to
% put large files on github.
% 
%    dataDirGitRepo = syncDataFromServer(projectPath)
%
% INPUTS
%   [projectPath]   :   project path on server (string)
%
% OUTPUTS
%   dataDirGitRepo  :   directory where data are synced to (string) 
%
%
% Example:
%   dataDirGitRepo = syncDataFromServer();


% Define where data live
if ~exist('projectPath', 'var') || isempty(projectPath)
    projectPath    = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
end

% Path to data on the server
dataDirServer  = fullfile(projectPath, 'data');

% Path to data on local machine (git repo)
dataDirGitRepo = fullfile(pfRV1rootPath, 'external');

% Make data directory if it does not exist
if ~exist(dataDirGitRepo, 'dir'); mkdir(dataDirGitRepo); end

% Sync data
system(sprintf('rsync -zvhr %s %s', dataDirServer, dataDirGitRepo));


return