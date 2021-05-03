function makeFigure7_Sensitivity_and_Asymmetry()
% Function to make Figure 7 of the manuscript:
% Radial asymmetries around the visual field: From retina to cortex to 
%   behavior. By Kupers, Benson, Carrasco, Winawer.
%    YEAR. JOURNAL. DOI.

%% 0. Define params and folders
nrSimConeDensities = 13;
c2rgc       = 0.5*(1:5).^2; % Cone 2 RGC ratios
expName     = 'conedensity';
meanPoissonPaddingFlag = true;
stimTemplateFlag       = false;

% Change folder names if using mean Poisson padded cone data
if meanPoissonPaddingFlag
    extraSubFolder = 'meanPoissonPadded';
else
    extraSubFolder = 'noPaddingBeforeConvolution';
end

if stimTemplateFlag
    subFolder = 'SVM-Energy';
    yl =  [1.3, 1.7];
    yl_asym = [-12, 50];
else
    subFolder = 'SVM-Fourier';
    yl =  [1.3, 1.7];
    yl_asym = [-10, 50];
end

% Folders
saveFigs     = true;
baseFolder   = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model';
dataFolder   = fullfile(baseFolder,'data',expName,'thresholds','rgc',extraSubFolder, subFolder,'average');
figureFolder = fullfile(baseFolder,'figures','surface3D_fixedRatios_rgc2c', expName, 'average', extraSubFolder,subFolder,'averge_lowess');
if (saveFigs) && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

%% 1. Load thresholds for each cone2mrgc ratio and fit function to data

% preallocate space
allData   = NaN(length(c2rgc), nrSimConeDensities);
errThresh = allData;
fct       = allData;
R2        = NaN(length(c2rgc),1);

for r = 1:length(c2rgc)
    
    % Load simulated contrast thresholds
    load(fullfile(dataFolder, sprintf('cThresholds_ratio%d_average.mat', r)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');
    allData(r,:) = [reshape(cell2mat(fit.ctrthresh),[],length(fit.ctrthresh))].*100;
    coneDensities  = xThresh;
    
    clear fit; clear xThresh fitToPlot dataToPlot;
    
    % Make "dummy" 2D grid
    [X,Y] = meshgrid(1:2,log10(coneDensities));
    Z = repmat(log10(allData(r,:))',1,2);
    idx = isfinite(Z);     % Find NaNs
    [meshFit, gof] = fit([X(idx) Y(idx)], Z(idx), 'lowess','span',0.3);
    
    % Extract single lines for separate ratio's
    fct(r,:) = meshFit(X(:,1),Y(:,1))';
    R2(r,:) = gof.rsquare;
end
 
%% Now fit a mesh to all points

% Find NaNs
idx = isfinite(allData');

% Make 2D grid
[X,Y] = meshgrid(log10(c2rgc),log10(coneDensities));
Z = log10(allData');

% Fit it!
[meshFit, gof] = fit([X(idx) Y(idx)], Z(idx), 'lowess','span',0.4);


%% 3. Predict contrast thresholds for given mRGC density
% if data does not exist, try syncDataFromServer()
% mRGC data for different meridia. Order=nasal, superior, temporal,inferior.
watson2014 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2');
watson2014.eccDeg = (0:0.05:40);

% Curcio et al. 1990 (left eye, retina coords, same order as mRGC)
curcio1990 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'conesCurcioISETBIO.mat'),'conesCurcioIsetbio', 'eccDeg','angDeg');
[~,angIdx] = intersect(curcio1990.angDeg,[0,90,180,270]);
coneDensityDeg2PerMeridian= curcio1990.conesCurcioIsetbio(angIdx,:);

% Compute RGC:cone ratio
rgc2coneRatio = watson2014.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2014.eccDeg==eccToCompute); % index
ratioAtIdx   = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Get cone density at chosen eccentricity for each meridian
observedConeDensitiesAtEccen = watson2014.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;

% Check: should be equal to curcio data
isequal(observedConeDensitiesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% take reciprocal for plotting -- meshfit expects cone2rgc ratio
ratioAtIdx = (1./ratioAtIdx);

% Find contrast threshold data for all meridians: Nasal, Superior,temporal, inferior
predictedContrastThreshold = 10.^meshFit(log10(ratioAtIdx),log10(observedConeDensitiesAtEccen));
predictedContrastThreshold = predictedContrastThreshold./100;

%% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
averageConeDensity_stimeccen = mean(observedConeDensitiesAtEccen);

for ii = 1:4
    errorRatioConeDensity(ii) = 2*abs(diff([observedConeDensitiesAtEccen(ii),averageConeDensity_stimeccen]));
end

% Get predicted thresholds for upper/lower error margins, using mesh fit
predictedErrorRGC = NaN(4,2);

% Nasal, superior, temporal, inferior retina (bounds: upper=1, lower=2) 
predictedErrorRGC(:,1) = 10.^meshFit(log10(ratioAtIdx),log10(observedConeDensitiesAtEccen'-errorRatioConeDensity));
predictedErrorRGC(:,2) = 10.^meshFit(log10(ratioAtIdx),log10(observedConeDensitiesAtEccen'+errorRatioConeDensity));
predictedErrorRGC = predictedErrorRGC./100;

%% Convert mean and error mRGC thresholds to sensitivity
rgcPredContrastSensitivityMEAN_retina        = 1./predictedContrastThreshold;
rgcPredContrastSensitivityERROR_retina_upperBound = 1./predictedErrorRGC(:,1);
rgcPredContrastSensitivityERROR_retina_lowerBound = 1./predictedErrorRGC(:,2);

retina2ToVisualFieldWithMeanHorz = @(x) [mean([x(1),x(3)]),x(4),x(2)];
visualField2Retina = @(x) [x(1),x(4),x(3),x(2)];
visualField2visualFieldWithMeanHorz = @(x) [mean([x(1),x(3)]),x(2),x(4)];

% Convert retinal coords into visual coords: HVM, UVM (inferior), LVM (superior)
rgcPredContrastSensitivityMEAN_wHorz_VF = retina2ToVisualFieldWithMeanHorz(rgcPredContrastSensitivityMEAN_retina);
rgcPredContrastSensitivityERROR_wHorz_VF_lower = retina2ToVisualFieldWithMeanHorz(rgcPredContrastSensitivityERROR_retina_lowerBound);
rgcPredContrastSensitivityERROR_wHorz_VF_upper = retina2ToVisualFieldWithMeanHorz(rgcPredContrastSensitivityERROR_retina_upperBound);

%% Get observed behavior + error
% OBSERVED Left HM, UVM, Right HM, LVM
obsContrastSensitivityMEAN_VF        = [46.4938; 28.9764; 47.7887; 34.3813]; % contrast senstivity (%)
obsContrastSensitivityERROR_VF       = [2.66468; 1.6445; 1.8450; 2.0505];    % contrast senstivity (%)
obsContrastSensitivityMEAN_retina    = visualField2Retina(obsContrastSensitivityMEAN_VF); %  L/R, LVM == superior retina, L/R, UVM == inferior retina

obsContrastSensitivityMEAN_wHorz_VF  = visualField2visualFieldWithMeanHorz(obsContrastSensitivityMEAN_VF);
obsContrastSensitivityERROR_wHorz_VF = visualField2visualFieldWithMeanHorz(obsContrastSensitivityERROR_VF);

%% Get simulated thresholds and error for cone model only 
% Cone data are from JWLOrientedGabor toolbox
% To get these data, run plotPsychometricFunctions('conedensity')
dataFolderConesOnly   = fullfile(baseFolder,'data',expName,'thresholds','absorptionsOnly',subFolder);
load(fullfile(dataFolderConesOnly,'coneabsorptionsOnly_predictedMeanAndError_stimeccen_linear'),'modelPredictionForPF','predictedError')
predictedErrorCones = predictedError; clear predictedError;

% MODELED CONE ONLY
conesPredContrastThreshMEAN_retina        = modelPredictionForPF; % nasal, superior, temporal, inferior
conesPredContrastThreshERROR_retina_lowerBound = predictedErrorCones(:,1); % - doubling diff in cone density from the mean
conesPredContrastThreshERROR_retina_upperBound = predictedErrorCones(:,2); % + doubling diff in cone density from the mean

% Get sensivity for retinal coords
conesPredContrastSensitivityMEAN_retina          = 1./conesPredContrastThreshMEAN_retina;

% Convert retinal coords to visual field coords
conesPredContrastSensitivityMEAN_wHorz_VF   = retina2ToVisualFieldWithMeanHorz(conesPredContrastSensitivityMEAN_retina);

% Do the same for upper/lower bounds of error margins
conesPredContrastSensitivityERROR_wHorz_VF_lowerBound = 1./retina2ToVisualFieldWithMeanHorz(conesPredContrastThreshERROR_retina_lowerBound);
conesPredContrastSensitivityERROR_wHorz_VF_upperBound = 1./retina2ToVisualFieldWithMeanHorz(conesPredContrastThreshERROR_retina_upperBound);
conesPredContrastSensitivityERROR_retina_lowerBound  = 1./conesPredContrastThreshERROR_retina_lowerBound; % - doubling diff in cone density from the mean
conesPredContrastSensitivityERROR_retina_upperBound   = 1./conesPredContrastThreshERROR_retina_upperBound; % + doubling diff in cone density from the mean

%% HVA VMA calc
HVAmean.obs       = hva(obsContrastSensitivityMEAN_retina);
VMAmean.obs       = vma(obsContrastSensitivityMEAN_retina);
HVAmean.predRGC   = hva(rgcPredContrastSensitivityMEAN_retina);
VMAmean.predRGC   = vma(rgcPredContrastSensitivityMEAN_retina);
HVAmean.predCones = hva(conesPredContrastSensitivityMEAN_retina);
VMAmean.predCones = vma(conesPredContrastSensitivityMEAN_retina);

HVAerror.obs       = HVAmean.obs + [-6.90, 6.90]; %  From Himmelberg et al. (2020) 
VMAerror.obs       = VMAmean.obs + [-5.65,5.65];  %  From Himmelberg et al. (2020)
HVAerror.predRGC   = [hva(rgcPredContrastSensitivityERROR_retina_lowerBound), hva(rgcPredContrastSensitivityERROR_retina_upperBound)];
VMAerror.predRGC   = [vma(rgcPredContrastSensitivityERROR_retina_lowerBound), vma(rgcPredContrastSensitivityERROR_retina_upperBound)];
HVAerror.predCones = [hva(conesPredContrastSensitivityERROR_retina_lowerBound), hva(conesPredContrastSensitivityERROR_retina_upperBound)];
VMAerror.predCones = [vma(conesPredContrastSensitivityERROR_retina_lowerBound), vma(conesPredContrastSensitivityERROR_retina_upperBound)];

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
errorbar(1:3,conesPredContrastSensitivityMEAN_wHorz_VF,...
             (conesPredContrastSensitivityMEAN_wHorz_VF-conesPredContrastSensitivityERROR_wHorz_VF_lowerBound), ...
             (conesPredContrastSensitivityERROR_wHorz_VF_upperBound-conesPredContrastSensitivityMEAN_wHorz_VF),...
             '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction Cones');

% Plot prediction for mRGCs
subplot(142)
bar(1:3, rgcPredContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,rgcPredContrastSensitivityMEAN_wHorz_VF,...
             (rgcPredContrastSensitivityMEAN_wHorz_VF-rgcPredContrastSensitivityERROR_wHorz_VF_lower), ...
            (rgcPredContrastSensitivityERROR_wHorz_VF_upper-rgcPredContrastSensitivityMEAN_wHorz_VF), ...
            '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction mRGCs');

% Plot prediction for behavior from Himmelberg, Winawer, Carrasco 2020 JoV
subplot(143)
bar(1:3, obsContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(3,:)); hold on
errorbar(1:3,obsContrastSensitivityMEAN_wHorz_VF,obsContrastSensitivityERROR_wHorz_VF, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
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
set(gca,'Xlim',[0,4],'Ylim',[yl_asym], 'TickDir', 'out', 'XTick', [1, 2.5], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Asymmetry');



if saveFigs
    hgexport(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    savefig(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    print(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'), '-dpng')
end


