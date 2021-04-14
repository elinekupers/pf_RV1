function [] = plotPsychometricFunctionsConeCurrent(baseFolder, expName, varargin)
% Function to compute psychometric functions based on the computational
% observer model.
%
% INPUTS:
% baseFolder         : path to project
% expName            : (string) the condition you want to plot.
%                       (See load expParams for possible conditions)
% [subFolder]        : (string) sub folder you want to plot data from.
% [saveFig]          : (boolean) save figures or not
% [plotAvg]          : (boolean) plot average across experiments runs or not
% [inputType]        : (string) use cone absorptions or cone current data
% [stimTemplateFlag] : (boolean) use SVM-Fourier (false) or stimulus based 
%                       energy template (true). Default is false. 
% [fitTypeName]      : (string) fitting function used to summarize cone 
%                       density vs thresholds: 
%                       for linearfit in log-log use 'linear',
%                       for robust linear fit in log-log (i.e. detecting 
%                       outliers and downweighting those) use 'linear-robust',
%                       for 2nd degree polynomial use 'poly2'. Default is
%                       'linear'.
%
% OUTPUTS:
% none
%
% Examples:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity', 'subFolder','run1')
% 
% plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity','subFolder','average','plotAvg',true)
%
% Using robust linear fit:
% plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity','subFolder','average','plotAvg',true, 'fitTypeName', 'robust-linear')
%% 0. Set general experiment parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addParameter('subFolder', @ischar);
p.addParameter('saveFig', true, @islogical);
p.addParameter('plotAvg', false, @islogical);
p.addParameter('inputType', 'current', @ischar);
p.addParameter('stimTemplateFlag',false, @islogical);
p.addParameter('fitTypeName','linear',@(x) ismember(x,{'linear','linear-robust','poly2'}));
p.parse(baseFolder, expName, varargin{:});

% Rename variables
baseFolder    = p.Results.baseFolder;
expName       = p.Results.expName;
subFolder     = p.Results.subFolder;
saveFig       = p.Results.saveFig;
plotAvg       = p.Results.plotAvg;
inputType     = p.Results.inputType;
stimTemplateFlag = p.Results.stimTemplateFlag;
fitTypeName   = p.Results.fitTypeName;

% Load specific experiment parameters
expParams    = loadExpParams(expName, false);
[xUnits, colors, labels, xThresh, ~] = loadWeibullPlottingParams('conedensity');

if stimTemplateFlag
    preFix = 'SVM-Energy';
else
    preFix = 'SVM-Fourier';
end

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification', inputType, preFix,subFolder);
figurePth   = fullfile(baseFolder,'figures','psychometricCurves', expName, inputType, [subFolder '_' preFix '_' fitTypeName]);
if ~exist(figurePth, 'dir')
    mkdir(figurePth); 
end

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = expParams.nTrials;

%% 1. Set Weibull fitting parameters

% Prepare fit variables
fit = [];

% Predefine cell arrays
fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

% Set inital slope, threshold for first stage fitting
fit.init   = [2, 0.01]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Get eccen at 4.5 deg;
% eccen  = 4.5;

% Check if different fitting accuracy for different cone types is requested

count = 1;
eccentricities = expParams.eccentricities;
for eccen = 1:length(eccentricities)
    
    %% 2. Get correct filename
    if plotAvg
        fName = sprintf('current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-0.60.30.1_AVERAGE.mat', eccentricities(eccen));
        fNameSE  = sprintf('current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-0.60.30.1_SE.mat', eccentricities(eccen));
        
        SE{count} = load(fullfile(dataPth, fNameSE));
        
    else
        
        if strcmp(expName, 'conedensity')
            fName = sprintf('current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00.mat', eccentricities(eccen));
        elseif strcmp(expName, 'defaultnophaseshift')
            fName = sprintf('current_Classify_coneOutputs_contrast0.1000_pa0_eye00_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat', eccentricities(eccen));
        end
    end

    %% 3. Load performance results
    
    % load model performance
    accuracy = load(fullfile(dataPth, fName));
    fn = fieldnames(accuracy);
    accuracy.P = squeeze(accuracy.(fn{1}));
    
    % Transpose matrix if necessary
    if size(accuracy.P,1)<size(accuracy.P,2)
        accuracy.P = accuracy.P';
    end
    
    if eccen == 1
        contrasts = expParams.contrastLevels;
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
    else
        contrasts = expParams.contrastLevelsPC;
        xUnits = linspace(min(expParams.contrastLevelsPC),max(expParams.contrastLevelsPC), 800);
    end
    
    %% 4. Fit Weibull
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, contrasts, accuracy.P(1:length(contrasts)), nTotal), fit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, xUnits);
    
    %% 5. Find contrast threshold
    fit.ctrthresh{count} = fit.ctrvar{count}(2);
    fit.data{count} = accuracy.P;
    
    count = count +1;
    
end % eccen


%% 6. Visualize psychometric curves

figure(3); clf; set(gcf,'Color','w', 'Position',  [218   316   986   488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s current', expName)); hold all;

% Remove empty cells
idx = ~cellfun(@isempty, fit.ctrpred);
fit.ctrpred = fit.ctrpred(idx);

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
if strcmp(expName, 'conedensity')
    logzero = 3e-3;
elseif strcmp(expName, 'defaultnophaseshift')
    logzero = 3e-5;
end

plotIdx = 1:length(fit.ctrpred);


% Loop over all functions to plot
for ii = plotIdx
    
    if ii == 1
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
    else
        xUnits = linspace(min(expParams.contrastLevelsPC),max(expParams.contrastLevelsPC), 800);
    end
    
    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;
    
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', '-');
    scatter(expParams.contrastLevelsPC(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
    
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
    if plotAvg
        errorbar([logzero, expParams.contrastLevelsPC(2:end)], dataToPlot, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevelsPC)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, 0.005, 0.01, 0.1, 1], 'XTickLabel',sprintfc('%1.1f',[0, 0.005, 0.01, 0.1, 1]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s_%s',expName,inputType, subFolder)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s_%s.eps',expName,inputType, subFolder)))
    print(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s_%s',expName,inputType, subFolder)), '-dpng')
end

if plotAvg
    % Save thresholds and fits
    thresholdsDir = fullfile(baseFolder,'data',expName,'thresholds', 'currentNoRGC',preFix);
    if ~exist(thresholdsDir, 'dir'); mkdir(thresholdsDir); end
    save(fullfile(thresholdsDir, sprintf('cThresholds_%s', subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh'); 


    % Plot density vs thresholds
    fNameSEThresh = sprintf('varThresh_coneResponse_current_13_conedensity.mat');
    load(fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC',preFix,fNameSEThresh),'varThresh');

    plotConeDensityVSThreshold(expName, fit, xThresh, 'varThresh', varThresh','inputType',inputType, 'fitTypeName',fitTypeName, 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log');    
end


%% Define error margins in terms of cone density, i.e.:
%% Get simulated thresholds and error for cone model only 
% Cone data are from JWLOrientedGabor toolbox
% To get these data, run plotPsychometricFunctions('conedensity')
source = fullfile(ogRootPath,'data','current',expName, sprintf('conecurrentOnly_predictedMeanAndError_stimeccen_%s.mat',fitTypeName));
destination = fullfile(thresholdsDir,sprintf('conecurrentOnly_predictedMeanAndError_stimeccen_%s.mat', fitTypeName));
eval(sprintf('!mv %s %s', source, destination));

load(fullfile(destination),'modelPredictionForPF','predictedError')

% MODELED CONE Current 
conesPredContrastThreshMEAN_retina        = modelPredictionForPF; % nasal, superior, temporal, inferior
conesPredContrastThreshERROR_retina_lower = predictedError(:,1); %#ok<NODEF> % - doubling diff in cone density from the mean
conesPredContrastThreshERROR_retina_upper = predictedError(:,2); % + doubling diff in cone density from the mean

% Get sensivity for retinal coords 
conesPredContrastSensitivityMEAN_retina          = 1./conesPredContrastThreshMEAN_retina;
% Convert retinal coords to visual field coords (superior/LVF swap with inferior/UVF)
conesPredContrastSensitivityMEAN_VF              = [conesPredContrastSensitivityMEAN_retina(1),conesPredContrastSensitivityMEAN_retina(4),conesPredContrastSensitivityMEAN_retina(3),conesPredContrastSensitivityMEAN_retina(2)];
% Combine left/right visual field coords to one horizontal value
conesPredContrastSensitivityMEAN_wHorz_VF        = [mean(conesPredContrastSensitivityMEAN_VF([1 3])),conesPredContrastSensitivityMEAN_VF(2),conesPredContrastSensitivityMEAN_VF(4)];

% Do the same for upper/lower bounds of error margins
conesPredContrastSensitivityERROR_wHorz_VF_lower = 1./[mean(conesPredContrastThreshERROR_retina_lower([1 3])),conesPredContrastThreshERROR_retina_lower(4),conesPredContrastThreshERROR_retina_lower(2)];
conesPredContrastSensitivityERROR_wHorz_VF_upper = 1./[mean(conesPredContrastThreshERROR_retina_upper([1 3])),conesPredContrastThreshERROR_retina_upper(4),conesPredContrastThreshERROR_retina_upper(2)];
conesPredContrastSensitivityERROR_retina_lower   = 1./conesPredContrastThreshERROR_retina_lower; % - doubling diff in cone density from the mean
conesPredContrastSensitivityERROR_retina_upper   = 1./conesPredContrastThreshERROR_retina_upper; % + doubling diff in cone density from the mean

%% HVA VMA calc

HVAmean = hva(conesPredContrastSensitivityMEAN_retina);
VMAmean = vma(conesPredContrastSensitivityMEAN_retina);
HVAerror = [hva(conesPredContrastSensitivityERROR_retina_lower), hva(conesPredContrastSensitivityERROR_retina_upper)];
VMAerror = [vma(conesPredContrastSensitivityERROR_retina_lower), vma(conesPredContrastSensitivityERROR_retina_upper)];

fprintf('HVA predicted for cone current\t %1.2f\n', HVAmean)
fprintf('VMA predicted for cone current\t %1.2f \n', VMAmean)
%%

% Bar plot to compare against behavior
condNames   = {'HVM', 'UVM','LVM'};
condColor   = [64, 224, 209]/255;

fH4 = figure(4); set(fH4, 'position',[383, 245, 1129, 542], 'color', 'w'); clf; 
subplot(121); hold all;
% Plot prediction for cones
bar(1:3, conesPredContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor); hold on
errorbar(1:3,conesPredContrastSensitivityMEAN_wHorz_VF,conesPredContrastSensitivityMEAN_wHorz_VF-conesPredContrastSensitivityERROR_wHorz_VF_lower,conesPredContrastSensitivityERROR_wHorz_VF_upper-conesPredContrastSensitivityMEAN_wHorz_VF,'.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[0.1, 1.4], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity'); title({'Model Prediction up to','cone phototransduction'});

% Plot Asymmetries in percent
subplot(122); hold on; cla
bar([0.5, 1], [HVAmean VMAmean], 0.2, 'EdgeColor','none','facecolor',condColor(1,:)); hold on

errorbar(0.5, HVAmean, HVAmean-HVAerror(1), HVAerror(2)-HVAmean, '.','color', 'k', 'LineWidth',2);
errorbar(1, VMAmean, diff([VMAmean,VMAerror(1)]), diff([VMAmean,VMAerror(2)]), '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0,1.5],'Ylim',[-20, 20], 'TickDir', 'out', 'XTick', [0.5, 1], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Polar angle asymmetry');


if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s_%s',expName, subFolder, fitTypeName)))
    hgexport(gcf,fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s_%s.eps',expName, subFolder, fitTypeName)))
    print(fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s_%s',expName, subFolder, fitTypeName)), '-dpng')
end
