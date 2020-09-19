function makeFigure7_Sensitivity_and_Asymmetry()
%
% Funxtion to make Figure 7 of the manuscript:
% TITLE. AUTHORS. JOURNAL. DOI.


%% 0. Define params and folders
nrSimConeDensities = 13;
c2rgc       = 1:5; % Cone 2 RGC ratios
rgc2c       = 2./c2rgc; % RGC 2 Cone ratios
expName     = 'conedensity';
subFolder   = 'average';
colors      = parula(length(c2rgc)+1);
labels      = sprintfc('RGC:cone = %1.1f:1.0', rgc2c);

% Folders
saveFigs     = true;
baseFolder   = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
dataFolder   = fullfile(baseFolder,'data',expName,'thresholds');
figureFolder = fullfile(baseFolder,'figures','surface3D', expName, subFolder);
if (saveFigs) && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

%% 1. Load thresholds for each cone2mrgc ratio and fit function to data

% preallocate space
allData   = NaN(length(c2rgc), nrSimConeDensities);
errThresh = allData;
fct       = allData;
R2        = NaN(length(c2rgc),1);

for r = c2rgc
    
    % Load simulated contrast thresholds
    load(fullfile(dataFolder, sprintf('cThresholds_ratio%d_%s', r, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');
    allData(r,:) = [reshape(cell2mat(fit.ctrthresh),[],length(fit.ctrthresh))].*100;
    
    % Load variance in simulated contrast thresholds computed by
    % bootstrapping simulation iterations starting with different rng seeds
    load(fullfile(dataFolder, sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptionrate_13_conedensity', r)), 'varThresh');
    errThresh(r,:) = varThresh.*100;
    
    % Fit a linear function in log-log space (so a powerlaw) for each simulated mrgc2cone ratio
    contrastThresh = cell2mat(fit.ctrthresh).*100; 
    coneDensities  = xThresh;
    lm = fitlm(log10(coneDensities),log10(contrastThresh));
    
    % Get fitted contrast thresholds (fct)
    fct(r,:) = lm.Coefficients.Estimate(2).*log10(xThresh) + lm.Coefficients.Estimate(1);
    
    % Get R2 of fits
    R2(r,:) = lm.Rsquared.ordinary;
end

% clean up 
clear fit xThresh fitToPlot dataToPlot;

%% 2. Upsample ratios and corresponding thresholds

dt = 41; % amount of upsampled values
c2rgc_upsampled  = linspace(c2rgc(1),c2rgc(end),dt);
rgc2c_upsampled  = 2./c2rgc_upsampled;

% Resample to log space, to get rectangles in mesh
coneDensities_resampled = 10.^linspace(log10(coneDensities(1)),log10(coneDensities(end)),13);

[X1,Y1] = meshgrid(c2rgc,coneDensities);
[X2,Y2] = meshgrid(c2rgc_upsampled,coneDensities_resampled);

fct_upsampled2 = griddata(X1,Y1, 10.^fct', X2,Y2, 'linear');

%% 3. Predict contrast thresholds for given mRGC density

% mRGC data for different meridia. Order=nasal, superior, temporal,inferior.
watson2015 = load(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2');
watson2015.eccDeg = (0:0.05:60);

% Curcio et al. 1990 (left eye, retina coords, same order as mRGC)
curcio1990 = load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'conesCurcioIsetbio', 'eccDeg','angDeg');
[~,angIdx] = intersect(curcio1990.angDeg,[0,90,180,270]);
coneDensityDeg2PerMeridian= curcio1990.conesCurcioIsetbio(angIdx,:);

% Compute RGC:cone ratio
rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2015.eccDeg==eccToCompute); % index
ratioAtIdx = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Find the ratio of interest in model grid for nasal and inferior retina
[err_rn,rn] = min(abs(rgc2c_upsampled-ratioAtIdx(1))); % nasal retina
[err_rs,rs] = min(abs(rgc2c_upsampled-ratioAtIdx(2))); % superior retina
[err_rt,rt] = min(abs(rgc2c_upsampled-ratioAtIdx(3))); % temporal retina
[err_ri,ri] = min(abs(rgc2c_upsampled-ratioAtIdx(4))); % inferior retina

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% Fit new line to upsampled contrast threshold data for all meridians
fct_upsampled2 = fct_upsampled2';
lm_rn = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rn,:))); % nasal
lm_rs = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rs,:))); % superior
lm_rt = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rt,:))); % temporal
lm_ri = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(ri,:))); % inferior

% Predicted threshold
cThreshold = @(x, a_coeff, b_intcpt) (a_coeff* log10(x)) + b_intcpt;

% Predicted cone density
% predictedDensity = @(y, a_coeff, b_intcpt) (10.^((y-b_intcpt)./a_coeff));

% Get contrast levels from model at 4.5 deg eccen
% Nasal retina
predictedContrastThreshold(1) = cThreshold(observedConesAtEccen(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedContrastThreshold(2) = cThreshold(observedConesAtEccen(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedContrastThreshold(3) = cThreshold(observedConesAtEccen(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedContrastThreshold(4) = cThreshold(observedConesAtEccen(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));

%% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
averageConeDensity_stimeccen = mean(observedConesAtEccen);

for ii = 1:4
    errorRatioConeDensity(ii) = 2*diff([observedConesAtEccen(ii),averageConeDensity_stimeccen]);
end

% Get predicted thresholds for those error margins, using their model fits
predictedThresholdDoubling = [];

% Nasal retina
predictedThresholdDoubling(1,1) = cThreshold(observedConesAtEccen(1)-errorRatioConeDensity(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
predictedThresholdDoubling(1,2) = cThreshold(observedConesAtEccen(1)+errorRatioConeDensity(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedThresholdDoubling(2,1) = cThreshold(observedConesAtEccen(2)-errorRatioConeDensity(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
predictedThresholdDoubling(2,2) = cThreshold(observedConesAtEccen(2)+errorRatioConeDensity(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedThresholdDoubling(3,1) = cThreshold(observedConesAtEccen(3)-errorRatioConeDensity(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
predictedThresholdDoubling(3,2) = cThreshold(observedConesAtEccen(3)+errorRatioConeDensity(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedThresholdDoubling(4,1) = cThreshold(observedConesAtEccen(4)-errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));
predictedThresholdDoubling(4,2) = cThreshold(observedConesAtEccen(4)+errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));

% Convert mean and error mRGC thresholds to sensitivity
rgcPredContrastSensitivityMEAN_retina        = 100.*(1./(10.^predictedContrastThreshold));
rgcPredContrastSensitivityERROR_retina       = 1./((10.^predictedThresholdDoubling)./100);
rgcPredContrastSensitivityERROR_retina_upper = rgcPredContrastSensitivityERROR_retina(:,1);
rgcPredContrastSensitivityERROR_retina_lower = rgcPredContrastSensitivityERROR_retina(:,2);

% Convert retinal coords into visual coords: HVM, UVM (inferior), LVM (superior)
rgcPredContrastSensitivityMEAN_wHorz_VF        = [mean(rgcPredContrastSensitivityMEAN_retina([1 3])); rgcPredContrastSensitivityMEAN_retina(4); rgcPredContrastSensitivityMEAN_retina(2)];
rgcPredContrastSensitivityERROR_wHorz_VF_lower = [mean(rgcPredContrastSensitivityERROR_retina_lower([1 3])); rgcPredContrastSensitivityERROR_retina_lower(4); rgcPredContrastSensitivityERROR_retina_lower(2)];
rgcPredContrastSensitivityERROR_wHorz_VF_upper = [mean(rgcPredContrastSensitivityERROR_retina_upper([1 3])); rgcPredContrastSensitivityERROR_retina_upper(4); rgcPredContrastSensitivityERROR_retina_upper(2)];

%% Get observed behavior + error
% OBSERVED LeftVM, UVM, RightVM, LVM
obsContrastSensitivityMEAN_VF        = [46.4938; 28.9764; 47.7887; 34.3813]; % contrast senstivity (%)
obsContrastSensitivityERROR_VF       = [2.66468; 1.6445; 1.8450; 2.0505];    % contrast senstivity (%)
obsContrastSensitivityMEAN_retina    = [obsContrastSensitivityMEAN_VF(1),obsContrastSensitivityMEAN_VF(4),obsContrastSensitivityMEAN_VF(3),obsContrastSensitivityMEAN_VF(2)];

obsContrastSensitivityMEAN_wHorz_VF  = [mean(obsContrastSensitivityMEAN_VF([1 3])), obsContrastSensitivityMEAN_VF(2), obsContrastSensitivityMEAN_VF(4)];
obsContrastSensitivityERROR_wHorz_VF = [mean(obsContrastSensitivityERROR_VF([1 3])), obsContrastSensitivityERROR_VF(2), obsContrastSensitivityERROR_VF(4)];


%% Get simulated thresholds and error for cone model only 
% Cone data are from JWLOrientedGabor toolbox
% To get these data, run plotPsychometricFunctions('conedensity')
load(fullfile(dataFolder,'coneOnly_predictedMeanAndError_stimeccen'),'modelPredictionForPF','predictedError')

% MODELED CONE ONLY
conesPredContrastThreshMEAN_retina        = modelPredictionForPF; % nasal, superior, temporal, inferior
conesPredContrastThreshERROR_retina_lower = predictedError(:,1); % - doubling diff in cone density from the mean
conesPredContrastThreshERROR_retina_upper = predictedError(:,2); % + doubling diff in cone density from the mean

conesPredContrastSensitivityMEAN_retina          =  1./conesPredContrastThreshMEAN_retina;
conesPredContrastSensitivityMEAN_VF              =  [conesPredContrastSensitivityMEAN_retina(1),conesPredContrastSensitivityMEAN_retina(4),conesPredContrastSensitivityMEAN_retina(3),conesPredContrastSensitivityMEAN_retina(2)];
conesPredContrastSensitivityMEAN_wHorz_VF        = [mean(conesPredContrastSensitivityMEAN_VF([1 3])),conesPredContrastSensitivityMEAN_VF(2),conesPredContrastSensitivityMEAN_VF(4)];
conesPredContrastSensitivityERROR_wHorz_VF_lower = 1./[mean(conesPredContrastThreshERROR_retina_lower([1 3])),conesPredContrastThreshERROR_retina_lower(4),conesPredContrastThreshERROR_retina_lower(2)];
conesPredContrastSensitivityERROR_wHorz_VF_upper = 1./[mean(conesPredContrastThreshERROR_retina_upper([1 3])),conesPredContrastThreshERROR_retina_upper(4),conesPredContrastThreshERROR_retina_upper(2)];
conesPredContrastSensitivityERROR_retina_lower   = 1./conesPredContrastThreshERROR_retina_lower; % - doubling diff in cone density from the mean
conesPredContrastSensitivityERROR_retina_upper   = 1./conesPredContrastThreshERROR_retina_upper; % + doubling diff in cone density from the mean

%% HVA VMA calc
HVAmean.obs       = hva(obsContrastSensitivityMEAN_retina);
VMAmean.obs       = vma(obsContrastSensitivityMEAN_retina);
HVAmean.predRGC   = hva(rgcPredContrastSensitivityMEAN_retina);
VMAmean.predRGC   = vma(rgcPredContrastSensitivityMEAN_retina);
HVAmean.predCones = hva(conesPredContrastSensitivityMEAN_retina);
VMAmean.predCones = vma(conesPredContrastSensitivityMEAN_retina);

HVAerror.obs       = HVAmean.obs + [-6.90, 6.90]; %  From Himmelberg et al. (2020) 
VMAerror.obs       = VMAmean.obs + [-5.65,5.65];  %  From Himmelberg et al. (2020)
HVAerror.predRGC   = [hva(rgcPredContrastSensitivityERROR_retina_lower), hva(rgcPredContrastSensitivityERROR_retina_upper)];
VMAerror.predRGC   = [vma(rgcPredContrastSensitivityERROR_retina_lower), vma(rgcPredContrastSensitivityERROR_retina_upper)];
HVAerror.predCones = [hva(conesPredContrastSensitivityERROR_retina_lower), hva(conesPredContrastSensitivityERROR_retina_upper)];
VMAerror.predCones = [vma(conesPredContrastSensitivityERROR_retina_lower), vma(conesPredContrastSensitivityERROR_retina_upper)];

combHVA = [HVAmean.predCones, HVAmean.predRGC, HVAmean.obs];
combVMA = [VMAmean.predCones, VMAmean.predRGC, VMAmean.obs];
errorCombHVA = [HVAerror.predCones(1), HVAerror.predRGC(1), HVAerror.obs(1); HVAerror.predCones(2), HVAerror.predRGC(2), HVAerror.obs(2)];
errorCombVMA = [VMAerror.predCones(1), VMAerror.predRGC(1), VMAerror.obs(1); VMAerror.predCones(2), VMAerror.predRGC(2), VMAerror.obs(2)];


fprintf('HVA predicted Cones \t %1.2f - predicted mRGC \t %1.2f - \t observed %1.2f \n', HVAmean.predCones, HVAmean.predRGC, HVAmean.obs)
fprintf('VMA predicted Cones \t %1.2f - predicted mRGC \t %1.2f - \t observed %1.2f \n', VMAmean.predCones, VMAmean.predRGC, VMAmean.obs)
%%

% Bar plot to compare against behavior

condNames   = {'HVM', 'UVM','LVM'};
condColor   = [63, 121, 204; 228, 65, 69;150 123 182]/255;

fH4 = figure(4); set(fH4, 'position',[383, 245, 1129, 542], 'color', 'w'); clf; hold all;
% Plot prediction for cones
subplot(141)
bar(1:3, conesPredContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(1,:)); hold on
errorbar(1:3,conesPredContrastSensitivityMEAN_wHorz_VF,conesPredContrastSensitivityMEAN_wHorz_VF-conesPredContrastSensitivityERROR_wHorz_VF_lower,conesPredContrastSensitivityERROR_wHorz_VF_upper-conesPredContrastSensitivityMEAN_wHorz_VF,'.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[1.3, 1.7], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction Cones');

% Plot prediction for mRGCs
subplot(142)
bar(1:3, rgcPredContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,rgcPredContrastSensitivityMEAN_wHorz_VF,rgcPredContrastSensitivityMEAN_wHorz_VF-rgcPredContrastSensitivityERROR_wHorz_VF_lower,rgcPredContrastSensitivityERROR_wHorz_VF_upper-rgcPredContrastSensitivityMEAN_wHorz_VF,'.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[1.3, 1.7], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction mRGCs');

% Plot prediction for behavior from Himmelberg, Winawer, Carrasco 2020 JoV
subplot(143)
bar(1:3, obsContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(3,:)); hold on
errorbar(1:3,obsContrastSensitivityMEAN_wHorz_VF,obsContrastSensitivityERROR_wHorz_VF, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[1.3, 1.7], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Behavior');

% Plot Asymmetries in percent
subplot(144); hold on; cla
x_bar = [0.5, 1, 1.5; 2.5, 3, 3.5];
for ii = 1:3
    bar(x_bar(:,ii), [combHVA(ii); combVMA(ii)], 0.2, 'EdgeColor','none','facecolor',condColor(ii,:)); hold on
end
errorbar(x_bar(1,:), combHVA, combHVA-errorCombHVA(1,:), errorCombHVA(2,:)-combHVA, '.','color', 'k', 'LineWidth',2);
errorbar(x_bar(2,:), combVMA, combVMA-errorCombVMA(1,:), errorCombVMA(2,:)-combVMA, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0,4],'Ylim',[-6, 50], 'TickDir', 'out', 'XTick', [1, 2.5], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Asymmetry');



if saveFigs
    hgexport(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    savefig(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    print(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'), '-dpng')
end


