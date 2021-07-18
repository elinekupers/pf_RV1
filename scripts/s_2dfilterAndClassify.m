% TODO:
% 1. check accuracy of classifier on absorptions and current for the no-eye
% movement L-only cone mosaic
% 2. how does performance change between the two, and is the differnce
% plausible?
% 3. Decide whether to start with absorptions of current for implementing
% filter + downsampling + noise
% 4 Run it and see if we get a systematic decline accuracy with the amount
% of downsampling
%
% See s_1dfilterAndClassify

pth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/conecurrent/defaultnophaseshift/';
expname = 'contrast0.1000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat';
runnum = 1;

load(fullfile(pth, sprintf('run%d', runnum), sprintf('current_OGconeOutputs_%s', expname)));
load(fullfile(pth, sprintf('run%d', runnum), sprintf('OGconeOutputs_%s', expname)));


% load the classification data
classify_absorptions = load(fullfile(pth, 'onlyL', sprintf('run%d', runnum), sprintf('Classify_coneOutputs_%s', expname)));
classify_currents    = load(fullfile(pth, 'onlyL', sprintf('run%d', runnum), sprintf('current_Classify_coneOutputs_%s', expname)));


figure, 
plot(classify_absorptions.expParams.contrastLevelsPC, classify_absorptions.accuracy,'o-', ...
    classify_absorptions.expParams.contrastLevelsPC, classify_currents.accuracy, ...
    'o-', 'LineWidth', 3, 'MarkerSize', 12); set(gca, 'XScale', 'log');

%% filter the data

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = false; % print figures or not
rgcParams.saveFigs   = false;
rgcParams.expName    = 'defaultnophaseshift';
rgcParams.subFolder  = [];

% Define DoG Params
ratio = 1;
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;               % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = 'absorptions';          % are we dealing with cone absorptions or current?
rgcParams.cone2RGCRatio = ratio;           % linear ratio
rgcParams.seed       = runnum;
cone2RGCRatio        = ratio;

% run the filter on the cone absorptions
data = mean(absorptions(:,:,:,1:28,:), 4);
[rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(data, rgcParams, expParams);


% next we need to 
% - classify the filtered data
% - then add noise 
% - then subsample at a few rates and classify again for each subsampling
% then if it makes sense, run with multiple contrasts on multiple runs


%% visualize the cone currents over time for one trial
figure(1)
for ii = 1:28
    subplot(1,2,1)
    jj= ii+10;
    imagesc(squeeze(current(1,:,:,jj,1)), [-26 -22]); 
    axis square; title(jj); 

    subplot(1,2,2)
    imagesc(squeeze(absorptions(1,:,:,ii,1)), [200 240]); 
    axis square; title(ii); waitforbuttonpress();%   pause(0.1); 
end

%%
noisyAbsorptions = absorptions; clear absorptions;

rng(1);
expParams.deg2m = 0.3 * 0.001;
expParams.polarAngle = -pi/2;
[cMosaic, cparams] = getConeMosaic(cparams.eccentricity, expParams, sparams);
cMosaic.integrationTime = 0.002; % ms

% Update pigment width/height
propCovered = getBanks1991ConeCoverage(cparams.eccentricity); % proportion
cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered; % meters
cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered; % meters

sparams.oi = oiDefocus(0, false); % input is 0 Zernicke defocus coeff
sparams.gabor.contrast  = 1.000;  % Michelson, range = [0 1]
sparams.freqCPD         = 4.00; % Cycles/degree
sparams.phases          = [pi/2 3*pi/2];
OG                      = ogStimuli(sparams);

cMosaic.noiseFlag = 'none';
n = 4; % use the 4th stimulus to go with the stored emPaths
noiselessAbsorptions = cMosaic.compute(OG(n), 'emPaths', emPaths, 'seed', 1+n);

% check noiseless vs noisy absorptions
figure(1), clf;
for ii = 1:28
    imagesc([squeeze(noiselessAbsorptions(1,:,:,ii)) squeeze(noisyAbsorptions(1,:,:,ii, 4))]);
    title(ii); colorbar; waitforbuttonpress();
end

figure(2), clf;
subplot(131); 
imagesc(squeeze(mean(noiselessAbsorptions(1,:,:,:),4)));
colorbar; title('Noiseless absorptions'); axis image
subplot(132); 
imagesc(squeeze(mean(noisyAbsorptions(1,:,:,:,4),4)));
colorbar; title('Loaded noisy absorptions'); axis image
subplot(133); 
imagesc(squeeze(mean(noiselessAbsorptions(1,:,:,:),4))-squeeze(mean(noisyAbsorptions(1,:,:,:,4),4)));
colorbar; title('Difference (the noise)'); axis image


%%

nTrials = size(noisyAbsorptions,1)*size(noisyAbsorptions,5);
rows = size(noisyAbsorptions,2);
cols = size(noisyAbsorptions,3);
nTimePoints = size(noisyAbsorptions,4);

workerID    = [];
meanCur     = [];
LMSfilters  = [];
currentSeed = 1;

photocurrents = zeros(nTrials, rows, cols, ...
    nTimePoints, 'single');
for ii = 1:2%nTrials
    if (~isempty(workerID)) && (mod(ii, round(nTrials / 10)) == 0)
        displayProgress(workerID, sprintf('%s-current', ...
            workDescription), 0.5 + 0.5 * ii / nTrials);
    end
    % Put this trial of absorptions into the cone mosaic
    cMosaic.absorptions = squeeze(noiselessAbsorptions(ii, :, :, :));
    currentSeed = currentSeed  + 1;
    if ii == 1 && (isempty(meanCur) || isempty(LMSfilters))
        % On the first trial, compute the interpolated linear
        % filters and the mean current
        [LMSfilters, meanCur] = cMosaic.computeCurrent('seed', ...
            currentSeed);
    else
        % On the remaining trials, use the same interpolated
        % filters and mean current
        LMSfilters = cMosaic.computeCurrent('seed', currentSeed, ...
            'interpFilters', LMSfilters, 'meanCur', meanCur);
    end
    photocurrents(ii, :, :, :) = single(reshape(cMosaic.current, ...
        [1 rows cols nTimePoints]));
end

%% can we get photocurrent output from absorption input that matches the stored photocurrent data?