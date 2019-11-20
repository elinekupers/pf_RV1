function [conesSong, dataMeridians] = getConesSong()
% Get cone density from Song data between 1-6 degrees eccentricity.
% These data were used on Noah Benson's SFN 2019 poster.

% Locate file
dataFile = fullfile(pfRV1rootPath, 'external', 'data', 'song_cone_density.mat');

if ~exist('dataFile')
    syncDataFromServer();
end

% Load data
conesSong = load(dataFile);

% Assuming the delta angle for each meridian is 15 deg.
eccen = (conesSong.eccentricity>=1) & (conesSong.eccentricity<=6); % deg

% % flip nasal / temporal
% n = conesSong.nasal;
% t = conesSong.temporal;
% i = conesSong.inferior;
% s = conesSong.superior;
% 
% conesSong.nasal = t;
% conesSong.temporal = n;
% conesSong.superior = i;
% conesSong.inferior = s;

integralPA15.nasal    = trapz(conesSong.nasal(eccen));
integralPA15.superior = trapz(conesSong.superior(eccen));
integralPA15.temporal = trapz(conesSong.temporal(eccen));
integralPA15.inferior = trapz(conesSong.inferior(eccen));

dataMeridians = [integralPA15.nasal, integralPA15.superior, integralPA15.temporal, integralPA15.inferior];

% fprintf('Cone density, 1-6 deg eccen, 15 deg polar angle wedge')
% fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(dataMeridians))
% fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(dataMeridians))

end