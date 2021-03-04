function [] = plotPsychometricFunctionsConeCurrent(baseFolder, expName, varargin)
% Function to compute psychometric functions based on the computational
% observer model.
%
% INPUTS:
% baseFolder      : path to project
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% [subFolder]     : string defining the sub folder you want to plot from.
% [saveFig]       : boolean defining to save figures or not
% [plotAvg]       : boolean defining to plot average across experiments runs or not
% [inputType]     : string defining data: absorptions or cone current
%
% OUTPUTS:
% none
%
% Examples:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity', 'subFolder','run1')
% 
% plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity','plotAvg',true)
%% 0. Set general experiment parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addParameter('subFolder', @ischar);
p.addParameter('saveFig', true, @islogical);
p.addParameter('plotAvg', false, @islogical);
p.addParameter('inputType', 'current', @ischar);
p.parse(baseFolder, expName, varargin{:});

% Rename variables
baseFolder    = p.Results.baseFolder;
expName       = p.Results.expName;
subFolder     = p.Results.subFolder;
saveFig       = p.Results.saveFig;
plotAvg       = p.Results.plotAvg;
inputType     = p.Results.inputType;

% Load specific experiment parameters
expParams    = loadExpParams(expName, false);
[xUnits, colors, labels, xThresh, ~] = loadWeibullPlottingParams('current');

for ll = 1:length(labels)
    tmp = strsplit(labels{ll}, ' ');
    polarAngleLabels{ll} = tmp{1};
    
    switch polarAngleLabels{ll}
        case 'nasal'
            polarAngles(ll) = 0;
        case 'superior'
            polarAngles(ll) = rad2deg(pi/2);
        case 'temporal'
            polarAngles(ll) = rad2deg(pi);
        case 'inferior'
            polarAngles(ll) = rad2deg(3*pi/2);
    end      
end

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification', inputType, expName);
figurePth   = fullfile(baseFolder,'figures','psychometricCurves', expName, [inputType 'allTimePoints']);
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
eccen  = 4.5;

% Check if different fitting accuracy for different cone types is requested

count = 1;
for pa = 1:length(polarAngleLabels)
    
    %% 2. Get correct filename
    if plotAvg
        subFolderName = sprintf('average_4.5deg_%s_currentAllTimePoints',polarAngleLabels{pa});
        fName = sprintf('current_Classify_coneOutputs_contrast1.000_pa%d_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_AVERAGE.mat', polarAngles(pa), eccen);
        fNameSE  = sprintf('current_Classify_coneOutputs_contrast1.000_pa%d_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_SE.mat', polarAngles(pa), eccen);
        
        SE{count} = load(fullfile(dataPth, subFolderName, fNameSE));
        
    else
        
        if strcmp(expName, 'conedensity')
            subFolderName = sprintf('%s_4.5deg_%s_currentAllTimePoints',subFolder, polarAngleLabels{pa});
%         fName = sprintf('current_Classify_coneOutputs_contrast1.000_pa%d_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00.mat', polarAngles(pa), eccen);
            fName = sprintf('current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00.mat', eccen);
        elseif strcmp(expName, 'defaultnophaseshift')
            subFolderName = sprintf('%s_currentAllTimePoints',subFolder);
            fName = sprintf('current_Classify_coneOutputs_contrast0.1000_pa0_eye00_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat', eccen);
        end
    end

    %% 3. Load performance results
    
    % load model performance
    accuracy = load(fullfile(dataPth,subFolderName, fName));
    fn = fieldnames(accuracy);
    accuracy.P = squeeze(accuracy.(fn{1}));
    
    % Transpose matrix if necessary
    if size(accuracy.P,1)<size(accuracy.P,2)
        accuracy.P = accuracy.P';
    end
    
%     xUnits = linspace(min(expParams.contrastLevelsPC),max(expParams.contrastLevelsPC), 100);
    
    %% 4. Fit Weibull
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevelsPC, accuracy.P, nTotal), fit.init);
    
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
% Only plot first two, 50:50 and last 2 functions for conetypes mixed experiment
plotIdx = 1:length(fit.ctrpred);


% Loop over all functions to plot
for ii = plotIdx

    [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('current');

    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;
    
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', lineStyles{ii});
    scatter(expParams.contrastLevelsPC(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
    
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
    if plotAvg
        errorbar([logzero, expParams.contrastLevelsPC(2:end)], dataToPlot, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevelsPC)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, expParams.contrastLevelsPC(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevelsPC(2:2:end)]*100))

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
    thresholdsDir = fullfile(baseFolder,'data',expName,'thresholds', 'currentNoRGC_allTimePoints',subFolderName);
    if ~exist(thresholdsDir, 'dir'); mkdir(thresholdsDir); end
    save(fullfile(thresholdsDir, sprintf('cThresholds_%s', subFolderName)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh'); 


    % Plot density vs thresholds
    fNameSEThresh = sprintf('varThresh_coneResponse_current_5_conedensity.mat');
    load(fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC_allTimePoints',fNameSEThresh),'varThresh');

    plotConeDensityVSThreshold(expName, fit, xThresh, 'varThresh', varThresh','inputType',inputType, 'fitTypeName','linear', 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log');    
end
%% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
% dataFolder   = fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC');
% load(fullfile(dataFolder,'conecurrentOnly_predictedMeanAndError_stimeccen'),'modelPredictionForPF','predictedError')
% 
% % MODELED CONE ONLY
% currentPredContrastThreshMEAN_retina        = cell2mat(fit.ctrthresh); % nasal, superior, temporal, inferior
% currentPredContrastThreshERROR_retina_lower = predictedError(:,1); %#ok<NODEF> % - doubling diff in cone density from the mean
% currentPredContrastThreshERROR_retina_upper = predictedError(:,2); % + doubling diff in cone density from the mean
% 
% % Get sensivity for retinal coords
% currentPredContrastSensitivityMEAN_retina          = 1./currentPredContrastThreshMEAN_retina;
% % Convert retinal coords to visual field coords
% currentPredContrastSensitivityMEAN_VF              = [currentPredContrastSensitivityMEAN_retina(1),currentPredContrastSensitivityMEAN_retina(4),currentPredContrastSensitivityMEAN_retina(3),currentPredContrastSensitivityMEAN_retina(2)];
% % Combine left/right visual field coords to one horizontal value
% currentPredContrastSensitivityMEAN_wHorz_VF        = [mean(currentPredContrastSensitivityMEAN_VF([1 3])),currentPredContrastSensitivityMEAN_VF(2),currentPredContrastSensitivityMEAN_VF(4)];
% 
% % Do the same for upper/lower bounds of error margins
% currentPredContrastSensitivityERROR_wHorz_VF_lower = 1./[mean(currentPredContrastThreshERROR_retina_lower([1 3])),currentPredContrastThreshERROR_retina_lower(4),currentPredContrastThreshERROR_retina_lower(2)];
% currentPredContrastSensitivityERROR_wHorz_VF_upper = 1./[mean(currentPredContrastThreshERROR_retina_upper([1 3])),currentPredContrastThreshERROR_retina_upper(4),currentPredContrastThreshERROR_retina_upper(2)];
% currentPredContrastSensitivityERROR_retina_lower   = 1./currentPredContrastThreshERROR_retina_lower; % - doubling diff in cone density from the mean
% currentPredContrastSensitivityERROR_retina_upper   = 1./currentPredContrastThreshERROR_retina_upper; % + doubling diff in cone density from the mean


%% Plot sensitivity and asymmetry in percent

condNames   = {'HVM', 'UVM','LVM'};

% Change to visual field coordinates and average left/right
thresholdsRetina  = cell2mat(fit.ctrthresh);
thresholdsVF       = [thresholdsRetina(1), ... HM
                     thresholdsRetina(4), ... inferior retina == UVM
                     thresholdsRetina(3), ... HM 
                     thresholdsRetina(2)]; % superior retina == LVM
                       
thresholdsVF_Horz = [mean(thresholdsVF([1 3])), ... HM
                           thresholdsVF(2), ... UVM
                           thresholdsVF(4)]; % LVM

% Get sensitivity from 1/thresholds
sensitivityMean_VFHorz  = 1./thresholdsVF_Horz;
sensitivityMean_Retina  = 1./thresholdsRetina;

if plotAvg
    % Use variance across simulation iterations to get errorbars  
    errRetina                = varThresh;
    errVF                    = [errRetina(1), ... HM
                                errRetina(4), ... inferior retina == UVM
                                errRetina(3), ... HM
                                errRetina(2)]; ... superior retina == LVM
    errVF_Horz               = [mean(errVF([1 3])), ... HM
                                errVF(4), ... inferior retina == UVM
                                errVF(2)]; ... superior retina == LVM
    sensitivityVFHorzErrorUpper  = 1./(thresholdsVF_Horz-errVF_Horz);
    sensitivityVFHorzErrorLower  = 1./(thresholdsVF_Horz+errVF_Horz);

    sensitivityRetinaErrorUpper      = 1./(thresholdsRetina-errRetina');
    sensitivityRetinaErrorLower      = 1./(thresholdsRetina+errRetina');
end

% Bar plot to compare against behavior
fH3 = figure(3); clf; set(fH3, 'position',[383   356   410   431], 'color', 'w');  hold all;

subplot(121)
bar(1:3, sensitivityMean_VFHorz,'EdgeColor','none','facecolor',colors(3,:)); hold on
if plotAvg
    errorbar(1:3,sensitivityMean_VFHorz,sensitivityMean_VFHorz-sensitivityVFHorzErrorLower,sensitivityVFHorzErrorUpper-sensitivityMean_VFHorz,'.','color', 'k', 'LineWidth',2);
end
set(gca,'Xlim',[0.2,3.8],'Ylim',[1, 20], 'TickDir', 'out', ...
    'XTick', [1:3], 'XTickLabel', condNames, ...
    'YTick', [1,10], 'YTickLabel', {'1', '10'}, ...
    'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model prediction cone current');

% Get HVA and VMA, input order is:
% (1) EAST (2) NORTH (3) WEST (4) SOUTH
HVAmean.current = hva(sensitivityMean_Retina);
VMAmean.current = vma(sensitivityMean_Retina);

if plotAvg
    HVAerrorUpper.current = hva(sensitivityRetinaErrorUpper);
    VMAerrorUpper.current = vma(sensitivityRetinaErrorUpper);
    HVAerrorLower.current = hva(sensitivityRetinaErrorLower);
    VMAerrorLower.current = vma(sensitivityRetinaErrorLower);
end

% Plot Asymmetries in percent
subplot(122); hold on; cla
x_bar = [0.5, 1];

bar(x_bar, [HVAmean.current; VMAmean.current], 0.2, 'EdgeColor','none','facecolor',colors(3,:)); hold on

if plotAvg
    errorbar(x_bar(1),HVAmean.current, HVAmean.current-HVAerrorLower.current, HVAerrorUpper.current-HVAmean.current, '.','color', 'k', 'LineWidth',2);
    errorbar(x_bar(2), VMAmean.current, VMAmean.current-VMAerrorLower.current, VMAerrorUpper.current-VMAmean.current, '.','color', 'k', 'LineWidth',2);
end

yl = minmax([HVAmean.current, VMAmean.current]);
yl(1) = yl(1)-5; yl(2) = yl(2)+5;

set(gca,'Xlim',[0,1.5],'Ylim',yl, 'TickDir', 'out', 'XTick', [0.5, 1], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Asymmetry');


if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s',expName, subFolder)))
    hgexport(gcf,fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s.eps',expName, subFolder)))
    print(fullfile(figurePth,sprintf('sensitivity_conecurrent_%s_%s',expName, subFolder)), '-dpng')
end
