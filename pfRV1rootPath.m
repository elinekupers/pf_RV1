function rootPath = pfRV1rootPath()
% Return the path to the root performance field Retina-V1 directory
%
% This function must reside in the directory at the base of the
% forward modeling code directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(pf_RV1_rootPath,'psychophys')

rootPath=which('pfRV1rootPath');

rootPath=fileparts(rootPath);

return
