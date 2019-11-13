function [rgcWatson, dataMeridians] = getMRGCRFWatson()
% Get midget RGC RF density from Watson data between 1-6 deg eccentricity.
% These data were used on Noah Benson's SFN 2019 poster.

%% Cones
rgcWatson = load(fullfile(pfRV1rootPath, 'external', 'data', 'watson_RGC_density.mat'));

% Assuming the delta angle for each meridian is 15 deg.
eccen = (rgcWatson.eccentricity>=1) & (rgcWatson.eccentricity<=6); % deg

integralPA15.nasal    = trapz(rgcWatson.nasal(eccen));
integralPA15.superior = trapz(rgcWatson.superior(eccen));
integralPA15.temporal = trapz(rgcWatson.temporal(eccen));
integralPA15.inferior = trapz(rgcWatson.inferior(eccen));

dataMeridians = [integralPA15.nasal, integralPA15.superior, integralPA15.temporal, integralPA15.inferior];

% fprintf('Midget RGC RF density, 1-6 deg eccen, 15 deg polar angle wedge')
% fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(dataIn))
% fprintf('Vertical-Meridian Asymmetry (North / South):\t %1.2f%%\n', vma(dataIn))

end