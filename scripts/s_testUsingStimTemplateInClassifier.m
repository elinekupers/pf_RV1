%% s_testUsingStimTemplateInClassifier.m
%
% Tinker script to test how adding a stimulus template before classifying
% RGC responses will affect SVM classifier accuracy.

% Try for defaultnophaseshift through RGC layer,
% i.e.: at 4.5 deg eccen, L-cones only, no eye movements, no phase
% differences within an oriented Gabor class, with Poisson noise at cone level.

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'defaultnophaseshift';
seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
ratio      = 1; % 2 = cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell)
eccen      = 1; % equal to 4.5 deg eccen

[PEnergy, PLinear] = linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen);

%% Visualize differences in psychometric functions

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'defaultnophaseshift';
eccen      = 1; % equal to 4.5 deg eccen

figure; set(gcf, 'Position',[ 332    39   618   759], 'color','w');
for ratio = 1:5
    % Load WITH stim template classifier results
    a=load(fullfile(baseFolder, 'data', expName, 'classification/rgc',  'stimTemplate',subFolder, ...
        sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_defaultnophaseshift_%s.mat',ratio,eccen,subFolder)));
    
    % Load WITHOUT stim template classifier results (note, includes all 5
    % ratio's in one file)
    b=load(fullfile(baseFolder, 'data', expName, 'classification/rgc', 'onlyL', ...
        'classifySVM_rgcResponse_Cones2RGC5_absorptionrate.mat'));
    
    contrastLevels.withTemplate = a.expParams.contrastLevels;
    contrastLevels.withoutTemplate = b.expParams.contrastLevels;
    
    performance.templateEnergy = a.P_svmEnergy;
    performance.templateLinear = a.P_svmLinear;
    performance.noTemplate = b.P_svm(ratio,:);
    
    logzero = 5.^-4;
    
    subplot(3,1,1)
    plot(contrastLevels.withTemplate, performance.templateEnergy, 'o-', 'lineWidth',2); hold on;
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    plot(logzero, data.withTemplate(1), 'o');
    title({sprintf('SVM-Energy Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})
    
    subplot(3,1,3)
    plot(contrastLevels.withTemplate, performance.templateLinear, 'o-', 'lineWidth',2); hold on;
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    plot(logzero, data.withTemplate(1), 'o');
    title({sprintf('SVM-Linear Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})
    
    
    subplot(3,1,3)
    plot(contrastLevels.withoutTemplate, data.withoutTemplate, 'o-', 'lineWidth',2); hold on;
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    plot(logzero, data.withoutTemplate(1), 'o')
    
    title({sprintf('SVM-Fourier Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})
 
    
end

% Add labels
for r = 1:5; labels{r} = sprintf('%dRGC2Cones',r); end
subplot(211);
l = findobj(gca);
legend(l([3:2:end]), labels, 'Location', 'NorthWest')
legend boxoff

subplot(212)
l = findobj(gca);
legend(l([3:2:end]), labels, 'Location', 'NorthWest')
legend boxoff

%% Try for cone density RGC simulations:
% Eccen: 0-40 deg, LMS cone mosaics, drift and microsaccades, 2 phases
% within an oriented Gabor class, and Poisson noise at cone level.

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'conedensity';
seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
ratio      = 1; % cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell)
eccen      = 5; % equal to 4.5 deg eccen

linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen)

%% Visualize differences in psychometric functions

figure; set(gcf, 'Position',[ 332    39   618   759], 'color','w');
for ratio = 1:5
    

    % Load WITH stim template classifier results
    a=load(fullfile(baseFolder, 'data', expName, 'classification','rgc', 'meanPoissonPadded', 'stimTemplate',subFolder, ...
        sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_%s.mat',ratio,eccen,subFolder)));
    
    % Load WITHOUT stim template classifier results
    b=load(fullfile(baseFolder, 'data', expName, 'classification/rgc', 'withPaddingBeforeConvolution', ...
        [subFolder '_meanPoissonPadded'], ...
        sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_%s_meanPoissonPadded.mat',ratio, eccen,subFolder)));
    
    
    contrastLevels.withTemplate = a.expParams.contrastLevels;
    contrastLevels.withoutTemplate = b.expParams.contrastLevels;
    
    data.withTemplate = a.P_svm;
    data.withoutTemplate = b.P_svm;
    
    subplot(2,1,1)
    plot(contrastLevels.withTemplate, data.withTemplate, 'o-'); hold on;
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    plot(10.^-3, data.withTemplate(1), 'o');
    title({sprintf('SVM Classifier performance for RGC responses WITH template'),'4.5 deg eccen, LMS-cones, eye movements, stimulus phase shifts'})
    
    
    subplot(2,1,2)
    plot(contrastLevels.withoutTemplate, data.withoutTemplate, 'o-'); hold on;
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    plot(10.^-3, data.withoutTemplate(1), 'o')
    
    title({sprintf('SVM Classifier performance for RGC responses WITHOUT template'),'4.5 deg eccen, LMS-cones, eye movements, stimulus phase shifts'})
    
end

% Add labels
for r = 1:5; labels{r} = sprintf('%dRGC2Cones',r); end
subplot(211);
l = findobj(gca);
legend(l([3:2:end]), labels, 'Location', 'NorthWest')
legend boxoff

subplot(212)
l = findobj(gca);
legend(l([3:2:end]), labels, 'Location', 'NorthWest')
legend boxoff