function dataDir = pfmf_InitDataDir

dataDir = fullfile(pfRV1rootPath, 'psychophys', 'data');

if ~exist(dataDir, 'dir'), mkdir(dataDir); end


return