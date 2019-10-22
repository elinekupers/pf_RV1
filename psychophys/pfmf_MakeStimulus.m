function pfmf_Im  = pfmf_MakeStimulus(stimParams, display)
% Construct image sequence for one trial of performance field modelfest
%  contrast experiment
% im = pfmf_MakeStimulus(stimulus, display)

sizeInPixels = 512;

% amplitude
stimParams.contrast              = stimParams.stimRange;
stimParams.colorDir              = [1 1 1];
stimParams.present               = 1; % relevant for pres / abs experiments

% temporal
stimParams.cyclesPerSecond       = 0;   % drift rate (0 = stationary)
stimParams.duration              = 0.25;	% in seconds
stimParams.temporalEnvelopeShape = 'gaussian';   % 'gaussian';%ß % or 'raisedcos' or 'gaussian' (try spread=duration/6)
stimParams.temporalSpread        = 0.0625;         % used only for gaussian temporal envelope

% spatial
stimParams.testPosition          = 'L'; % irrelevant if ecc == 0 deg
stimParams.cyclesPerDegree       = 4;	% spatial frequency
stimParams.orientDegrees         = 0;  % horizontal
stimParams.phaseDegrees          = 0;
stimParams.size                  = 4.2667; % hard edge of stimulus in deg
stimParams.spread                = fwhm2sd(0.5); % gaussian sigma in deg


stimParams.eccentricity          = 4; % stimulus location in deg (negative = down)
stimParams.polarangle            = 0; % stimulus location in deg (negative = down)

% position of fixation cross
stimParams.fixationEcc           = 0; 

return

function sigma = fwhm2sd(fwhm)
% convert a gaussian fullwidth-halfmax to sd

sigma = fwhm / 2.35;

return