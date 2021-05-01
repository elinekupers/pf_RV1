function mRGCRFDensityPerDeg2 = getMRGCRFWatson(eccDeg)
%% Function to get midget RGC RF density from Watson 2014, along requested
% degrees of eccentricity. This function relies on the ISETBIO toolbox.

%% Instantiate a WatsonRGCModel object. 
% Set the 'generateAllFigures' flag to true to generate several figures of
% the Watson 2014 paper
WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);

% WARNING: Be aware that all figures and internal functions in the
% WatsonRGCmodel object assume coordinates as in the paper: visual field
% coordinates of the right eye. We, on the other hand, want meridians in
% retinal coordinates. So, now you might think. Hold on a second, shouldn't
% superior and inferior meridians be flipped then in your meridianLabel
% variable?! And then I'd say: Excellent observation! Fortunately, we don't
% have to use most of this model and can simply by-pass any assumed order
% of the meridians and flipping of axes by directing calling the function,
% mRGCRFSpacingAndDensityAlongMeridian. 
% If you don't believe me, go to line 216 in the WatsonRGCModel to check
% out the wrapper function called:
% mRGCRFSpacingAndDensityAtRetinalPositions. Here, you remap the assumed
% visual field meridian order to retinal coordinates before feeding it into
% the function that does the work mRGCRFSpacingAndDensityAlongMeridian.

% Compute total RGC RF density along all cardinal meridians for requested
% eccentricities.
meridianLabel = {'nasal meridian','superior meridian','temporal meridian','inferior meridian'};
for ii = 1:length(meridianLabel)
    [~,mRGCRFDensityPerDeg2(ii,:)] = WatsonRGCCalc.mRGCRFSpacingAndDensityAlongMeridian(eccDeg, meridianLabel{ii}, 'deg', 'deg^2');
end

end