function snr_prctCorrect = getSNRAccuracy(data, expParams, inputType, cone2RGCRatio, c, baseFolder)

% Check inputs
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC;
    selectTimePoints = 51:78;
end

% Find ideal observer (noiseless) RGC filterd cone data
idealObserverDir = fullfile(baseFolder, 'data', 'idealobserver', 'rgc', 'meanPoissonPadded', 'onlyL', sprintf('ratio%d',cone2RGCRatio));

% preallocate space
% snr_dB      = NaN(size(expParams.contrastLevels));
% snr_power    = NaN(size(expParams.contrastLevels));
% snr_prctCorrect = NaN(size(expParams.contrastLevels));


% for c = 1:length(expParams.contrastLevels)

idealObserverFile = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen4.50_%s.mat', cone2RGCRatio, contrasts(c), inputType);

templateNoiseless = load(fullfile(idealObserverDir, idealObserverFile));
theSignal = templateNoiseless.rgcResponse;
theNoise  = data - theSignal;

% Get the trials and samples (should be the data for all data sets though
nStimuli = size(theNoise,5);
nTrials  = size(theNoise,1) * nStimuli/2;
tSamples = size(theNoise,4);
nrows    = size(theNoise,2);
ncols    = size(theNoise,3);

% Label clockwise and counterclockwise trials
label = [ones(nTrials, 1); -ones(nTrials, 1)]; % First set is CW, second set is CCW

% Permute to trials x stimuli x rows x cols x time points
theNoise  = permute(theNoise, [1 5 2:4]);
theSignal = permute(theSignal, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
theNoise  = reshape(theNoise, [], nrows, ncols, tSamples);
theSignal = reshape(theSignal, [], nrows, ncols, tSamples);

% Define signal subtracting the two stimulus conditions
signalCW      = theSignal(label==1,:,:,selectTimePoints);
signalCCW     = theSignal(label==-1,:,:,selectTimePoints);
theSignalDiff = signalCW-signalCCW;
theSignalDiff = theSignalDiff(:);

% Get all trials with the same label and only take one trial and one time point
% (since data are all the same across time, without any photon noise or eyemovement or
% phase shifts)
theNoiseCW = theNoise(label==1,:,:,selectTimePoints);
theNoiseCW = theNoiseCW(:);

% Get SNR
snr_dB = snr(theSignalDiff, theNoiseCW);

% Convert SNR dB to dprime / power
snr_power = db2mag(snr_dB);

% Convert SNR dprime / power to percent correct
snr_prctCorrect = normcdf(snr_power/2);

fprintf('Contrast %1.4f \t SNR dB: %2.3f, d-prime/power: %2.3f,  percent correct: %2.3f\n', contrasts(c), snr_dB, snr_power, snr_prctCorrect)

return






