function rgcWatson = getMRGCRFWatson(eccDeg)
% Get midget RGC RF density from Watson data between 0-40 deg eccentricity.

%% %% -----------------------------------------------------------------
%  -------------- mRGC from Watson 2015 using ISETBIO --------------
%  -----------------------------------------------------------------

% Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag 
% to true to generate several figures of the Watson 2014 paper
WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);

% Compute total RGC RF density along the superior meridian for a number of eccentricities
meridianLabel = {'nasal meridian','superior meridian','temporal meridian','inferior meridian'};
for ii = 1:length(meridianLabel)
    rgcWatson(ii,:) = WatsonRGCCalc.midgetRGCRFDensity(eccDeg, meridianLabel{ii}, 'RFs per deg2');
end


end