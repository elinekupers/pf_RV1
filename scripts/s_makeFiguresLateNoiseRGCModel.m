%% s_plotFiguresLateNoiseRGCmodel

%% 0. Set general parameters
runnum            = 2;
pth               = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/';
expName           = 'conedensitynophaseshiftlonly500';
subfolder         = sprintf('run%d', runnum);
expParams         = loadExpParams(expName);
saveFig           = true;
lateNoiseLevel    = 1;
dataTypeLabels    = {'absorptions', 'current', 'Filtered', 'LateNoise',...
                     'DownSampled1', 'DownSampled2', 'DownSampled3',...
                     'DownSampled4', 'DownSampled5'};
downsampleFactors = 2./(1:5).^2; % RGC:cone downsample ratios for 2D arrays

% Convert eccentricity to cone density using Curcio et al. 1990 data
[conedensityLabels, cDensity] =  getConeDensityLabelsForPlotting(expParams);

if saveFig
    figurePth  = fullfile('/Volumes/server/Projects/PerformanceFields_RetinaV1Model/', ...
        'figures','psychometricCurves', expName, 'current', subfolder);
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
end

%% 1. Set Weibull parameters

nrDataTypes = length(dataTypeLabels);
nrEccen     = length(expParams.eccentricities);
colors      = jet(nrEccen+1);

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-4;

% Finely sample contrast levels to extract contrast threshold from Weibull
xUnits      = logspace(log10(logzero),log10(max(expParams.contrastLevels)), 500);

% Predefine arrays for Weibull fits and params
weibullFit = struct();
weibullFit.ctrvar       = NaN(nrDataTypes,nrEccen, 2); % estimated variables
weibullFit.ctrpred      = NaN(nrDataTypes,nrEccen, length(xUnits)); % estimated fine sampled Weibull prediction
weibullFit.data         = NaN(nrDataTypes,nrEccen, length(expParams.contrastLevels));
weibullFit.ctrthresh    = NaN(nrDataTypes,nrEccen, 1); % final threshold
weibullFit.init         = [4, 0.001]; % slope (at ~80%) and initial threshold

%% 2. Fit!
for eccen = 1:nrEccen
    loadStr = sprintf('classifierAccuracy_latenoiselevel%1.1f_eccen%1.2f_withDownsampling.mat', ...
        lateNoiseLevel, expParams.eccentricities(eccen));
    load(fullfile(pth,'conecurrent', expName, subfolder, 'classifierAccuracy',loadStr), 'PercentCorrect')
    
    % Loop over data types
    for dt = 1:nrDataTypes
        
        % Absorptions have a much lower threshold and steeper slope, 
        % so we change initial start values for better fit
        if dt == 1, weibullFit.init = [4, 0.001];
        else, weibullFit.init = [3, 0.01]; end
        
        accuracy = PercentCorrect(:,dt);
        contrasts = expParams.contrastLevels;

        % Find NaNs and remove them from input data 
        outliers = find(isnan(accuracy));
        if ~isempty(outliers)
            contrasts(outliers) = [];
            accuracy(outliers) = [];
        end
        
        % Fit weibul!
        w = fitWeibullToClassifierAccuracy(contrasts, accuracy, weibullFit.init, xUnits);
        
        % Store params
        weibullFit.ctrvar(dt,eccen, :) = w.ctrvar;
        weibullFit.ctrpred(dt, eccen, :) = w.ctrpred;
        weibullFit.data(dt, eccen,:) = PercentCorrect(:,dt);
        weibullFit.ctrthresh(dt,eccen,:) = w.ctrthresh;
    end
end

%% 3. Visualize weibulls per cone density and datatype

for dt = 1:nrDataTypes
    dataPoints = squeeze(weibullFit.data(dt, :,:));
    fittedWeibull = squeeze(weibullFit.ctrpred(dt,:,:));
    
    fH1(dt) = makeFigure_WeibullLateNoiseRGCModel(dataPoints, fittedWeibull, ...
         xUnits, expName, expParams, dataTypeLabels{dt}, saveFig, figurePth);
end

%% 4. Plot contrast thresholds vs downsampling

fH2 = makeFigure_ThreshVsDownsamplingLateNoiseRGCModel(...
    weibullFit, expName, expParams, figurePth, saveFig);

%% 5. Plot contrast thresholds vs cone density

whichfit = 'lowess';
[fH3, rgc3D, current3D, absorptions3D] = ...
    makeFigure_ThreshVsConeDensityLateNoiseRGCModel(...
    weibullFit, whichfit, dataTypeLabels, expName, expParams, figurePth, saveFig);

%% 6. Get cone density and down sample factor at 4.5 deg
% i.e, the eccentricity of the psychophysical experiment we compare model predictions to

% Get mRGC data for different meridia.
% Order = nasal, superior, temporal,inferior.
watson2015 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'mRGCWatsonISETBIO.mat'), ...
    'mRGCRFDensityPerDeg2', 'eccDeg');
assert([length(watson2015.eccDeg) == length(0:0.05:40)]);

% Curcio et al. 1990 (left eye, retina coords, same order as mRGC)
curcio1990 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'conesCurcioISETBIO.mat'), ...
    'conesCurcioIsetbio', 'eccDeg','angDeg');
[~,angIdx] = intersect(curcio1990.angDeg,[0,90,180,270]);
coneDensityDeg2PerMeridian= curcio1990.conesCurcioIsetbio(angIdx,:);

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2015.eccDeg==eccToCompute); % index

% Compute cone:RGC ratio
rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;
ratioAtIdx   = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% Find contrast threshold data for all meridians: Nasal, Superior,temporal, inferior
predictedContrastThreshold = 10.^rgc3D.meshFit(log10(ratioAtIdx),log10(observedConesAtEccen));

%% 7. Make 3D mesh plot
fH4 = makeFigure_3DSurfaceMesh_LateNoiseRGCModel(...
            rgc3D, ratioAtIdx, observedConesAtEccen, ...
            predictedContrastThreshold, expName, figurePth, saveFig);

%% 8. Get thresholds, sensitivity, error bars, asymmetries from fits

% Convert predicted thresholds from log10 fraction to fraction
% Retinal coords: % nasal, superior, temporal, inferior
prediction_retina.rgc.cThresholds.mean      = predictedContrastThreshold;
prediction_retina.cones.cThresholds.mean    = 10.^absorptions3D.meshFit(ones(4,1),log10(observedConesAtEccen));
prediction_retina.current.cThresholds.mean  = 10.^current3D.meshFit(ones(4,1),log10(observedConesAtEccen));

% Convert thresholds to sensitivity
prediction_retina.rgc.sensitivity.mean      = 1./prediction_retina.rgc.cThresholds.mean;
prediction_retina.cones.sensitivity.mean    = 1./prediction_retina.cones.cThresholds.mean;
prediction_retina.current.sensitivity.mean  = 1./prediction_retina.current.cThresholds.mean;

% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
averageConeDensity_stimeccen = mean(observedConesAtEccen);
diffConeDensityFromMean      = averageConeDensity_stimeccen-observedConesAtEccen;
errorRatioConeDensity        = 0.5*diffConeDensityFromMean;

% Get error for RGC responses, absorptions and current.
% Note: fits are in log10-log10, hence we convert inputs and outputs
prediction_retina.rgc.cThresholds.error.upper     = 10.^rgc3D.meshFit(log10(ratioAtIdx),log10(observedConesAtEccen+errorRatioConeDensity));
prediction_retina.rgc.cThresholds.error.lower     = 10.^rgc3D.meshFit(log10(ratioAtIdx),log10(observedConesAtEccen-errorRatioConeDensity));
prediction_retina.cones.cThresholds.error.upper   = 10.^absorptions3D.meshFit(ones(4,1),log10(observedConesAtEccen+errorRatioConeDensity));
prediction_retina.cones.cThresholds.error.lower   = 10.^absorptions3D.meshFit(ones(4,1),log10(observedConesAtEccen-errorRatioConeDensity));
prediction_retina.current.cThresholds.error.upper = 10.^current3D.meshFit(ones(4,1),log10(observedConesAtEccen+errorRatioConeDensity));
prediction_retina.current.cThresholds.error.lower = 10.^current3D.meshFit(ones(4,1),log10(observedConesAtEccen-errorRatioConeDensity));

% Convert retinal coords into visual coords:
% Retina to visual field, where we average nasal/retina:
%   HVM (average nasal & temporal), UVM (inferior), LVM (superior)
r2VF_wMeanHorz  = @(data) [mean([data(1),data(3)]),data(4),data(2)];

% Visual field to visual field where we average nasal/retina
%   left, upper, right, lower VF --> horizontal, upper, lower VF
vf2vf_wMeanHorz = @(data) [mean([data(1),data(3)]),data(2),data(4)];

% Convert error mRGC thresholds to sensitivity and visual field coords
modelStages = {'cones','current','rgc'};
for ii = 1:length(modelStages)
    prediction_retina.(modelStages{ii}).sensitivity.error.upper = 1./prediction_retina.(modelStages{ii}).cThresholds.error.upper;
    prediction_retina.(modelStages{ii}).sensitivity.error.lower = 1./prediction_retina.(modelStages{ii}).cThresholds.error.lower;
    
    prediction_visualfield.(modelStages{ii}).sensitivity.mean_wHorz = r2VF_wMeanHorz(prediction_retina.(modelStages{ii}).sensitivity.mean);
    prediction_visualfield.(modelStages{ii}).sensitivity.error.upper_wHorz = r2VF_wMeanHorz(prediction_retina.(modelStages{ii}).sensitivity.error.upper);
    prediction_visualfield.(modelStages{ii}).sensitivity.error.lower_wHorz = r2VF_wMeanHorz(prediction_retina.(modelStages{ii}).sensitivity.error.lower);
end

% Get observed behavior + error from psychophysics experiment
% OBSERVED HM, UVM, Right HM, LVM
% (data from baseline experiment Himmelberg, Winawer, Carrasco, 2020, JoV)
observed_visualfield.sensitivity.mean  = [46.4938; 28.9764; 47.7887; 34.3813]; % contrast senstivity (%)
observed_visualfield.sensitivity.error = [2.66468; 1.6445; 1.8450; 2.0505];    % contrast senstivity (%)

% Visual field to retina, so flip upper/lower VF, no averaging of horizontal
%   left, upper, right, lower VF --> nasal, superior, temporal, inferior
vf2r = @(data) [data(1),data(4),data(3),data(2)];
observed_retina.sensitivity.mean  = vf2r(observed_visualfield.sensitivity.mean); %  L/R, LVM == superior retina, L/R, UVM == inferior retina

observed_visualfield.sensitivity.mean_wHorz  = vf2vf_wMeanHorz(observed_visualfield.sensitivity.mean);
observed_visualfield.sensitivity.error_wHorz = vf2vf_wMeanHorz(observed_visualfield.sensitivity.error);

% Asymmetry (HVA VMA) calculations
HVAmean.obs         = hva(observed_retina.sensitivity.mean);
VMAmean.obs         = vma(observed_retina.sensitivity.mean);

HVAerror.obs       = HVAmean.obs + [-6.90, 6.90]; %  From Himmelberg et al. (2020)
VMAerror.obs       = VMAmean.obs + [-5.65,5.65];  %  From Himmelberg et al. (2020)

modelStages = {'cones','current','rgc'};
for ii = 1:length(modelStages)
    HVAmean.(modelStages{ii}) = hva(prediction_retina.(modelStages{ii}).sensitivity.mean);
    VMAmean.(modelStages{ii}) = vma(prediction_retina.(modelStages{ii}).sensitivity.mean);
    
    HVAerror.(modelStages{ii})   = [hva(prediction_retina.(modelStages{ii}).sensitivity.error.lower), ...
        hva(prediction_retina.(modelStages{ii}).sensitivity.error.upper)];
    VMAerror.(modelStages{ii})   = [vma(prediction_retina.(modelStages{ii}).sensitivity.error.lower), ...
        vma(prediction_retina.(modelStages{ii}).sensitivity.error.upper)];
    
    fprintf('HVA predicted for %s \t %1.2f%% [%1.2f%% %1.2f%%]\n', modelStages{ii}, ...
        HVAmean.(modelStages{ii}), HVAerror.(modelStages{ii}))
    fprintf('VMA predicted for %s \t %1.2f%% [%1.2f%% %1.2f%%]\n', modelStages{ii}, ...
        VMAmean.(modelStages{ii}), VMAerror.(modelStages{ii}))
    fprintf('\n')
end

combHVA = [HVAmean.cones, HVAmean.current, HVAmean.rgc, HVAmean.obs];
combVMA = [VMAmean.cones, VMAmean.current, VMAmean.rgc, VMAmean.obs];

errorCombHVA = [HVAerror.cones(1), HVAerror.current(1), HVAerror.rgc(1), HVAerror.obs(1); ...
    HVAerror.cones(2), HVAerror.current(2), HVAerror.rgc(2), HVAerror.obs(2)];

errorCombVMA = [VMAerror.cones(1), VMAerror.current(1), VMAerror.rgc(1), VMAerror.obs(1); ...
    VMAerror.cones(2), VMAerror.current(2), VMAerror.rgc(2), VMAerror.obs(2)];

%% 9. Plot sensitivity for each modeling stage + observed behavior as bars 
fH5 = makeFigure_PredictedSensivitiy_LateNoiseRGCModel(...
        prediction_visualfield, observed_visualfield, combHVA, errorCombHVA, combVMA, errorCombVMA, ...
        expName,figurePth,saveFig);


