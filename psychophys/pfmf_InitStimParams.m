function stimParams  = pfmf_InitStimParams(display, whichExperiment)

display.devices        = getDevices;
stimParams.inputDevice = getBestDevice(display);

% initialize stimulus parameters

% amplitude
stimParams.contrast              = 1.0;
stimParams.colorDir              = [1 1 1];
stimParams.present               = 1; % relevant for pres / abs experiments

% temporal
stimParams.cyclesPerSecond       = 0;   % drift rate (0 = stationary)
stimParams.duration              = 0.25;	% in seconds
stimParams.preStimduration       = 0.5;	% in seconds
stimParams.temporalEnvelopeShape = 'gaussian';   % 'gaussian';%ß % or 'raisedcos' or 'gaussian' (try spread=duration/6)
stimParams.temporalSpread        = 0.0625;       % used only for gaussian temporal envelope

% spatial
if strcmp(whichExperiment, 'East')
    stimParams.testPosition          = 0;    % irrelevant if ecc == 0 deg
elseif strcmp(whichExperiment == 'South') 
   stimParams.testPosition          = 90;    % irrelevant if ecc == 0 deg
elseif strcmp(whichExperiment == 'West')    
   stimParams.testPosition          = 180;   % irrelevant if ecc == 0 deg
elseif strcmp(whichExperiment == 'North') 
   stimParams.testPosition          = 270;  % irrelevant if ecc == 0 deg
end

stimParams.testPosition          = 0;   % [0 90 180 270]; % irrelevant if ecc == 0 deg
stimParams.cyclesPerDegree       = 4;	% spatial frequency
stimParams.orientDegrees         = 0;
stimParams.phaseDegrees          = 0;
stimParams.size                  = 4.26; % hard edge of stimulus in deg
stimParams.spread                = 0.5;% fwhm2sd(0.5); % gaussian sigma in deg

stimParams.eccentricity          = 5; % stimulus location in deg (negative = down)

% position of fixation cross
stimParams.fixationEcc           = 0; 


function sigma = fwhm2sd(fwhm)
% convert a gaussian fullwidth-halfmax to sd

sigma = fwhm / 2.35;

return



%% fixed parameters

% display.devices        = getDevices;
% stimParams.inputDevice = getBestDevice(display);
% 
% stimParams.exp         = whichExperiment;
% 
% % screen diameter in pixels
% screensize = min(display.numPixels);
% 
% % Stim params
% stimParams.duration  = 0.25; % seconds
% stimParams.cueduration   = 0.5; % seconds
% 
% % stimParam.nTimePoints = ???
% stimParams.isi      = 2;
% 
% % specify range of contrast values
% stimParams.contrastRange = [0 logspace(-4,-1.4,29)];
% stimParams.neutral   = zeros(length(stimParams.contrastRange));
% 
% % speficy range of contrast values
% % bk = display.backColorIndex;
% % mx = max(display.stimRgbRange);
% % tmp = display.backColorIndex + (1:(mx-bk));
% % contrast = tmp ./ display.backColorIndex - 1;
% % stimParams.con = contrast;
% 
% 
% if strcmp(whichExperiment, 'uniform')
%     stimParams.background  = []; 
% elseif strcmp(whichExperiment, '1/f')
%     stimParams.background  = [];
% end


%% debug


% % *** Pause Event parameters ***
% stimParams.fontSize = 18;
% stimParams.fontFace = 'Arial';
% stimParams.curStr = 'Take a break';
% stimParams.centerLoc = screensize/2;
% stimParams.angle = 0;

return