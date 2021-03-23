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

% [PEnergy, PLinear] = linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen);

%% Visualize differences in psychometric functions

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'defaultnophaseshift';
eccen      = 1; % equal to 4.5 deg eccen
cmap = parula(6);
figure; set(gcf, 'Position',[ 332    39   618   759], 'color','w');

% Load WITHOUT stim template classifier results (note, includes all 5
% ratio's in one file)
woTemplate=load(fullfile(baseFolder, 'data', expName, 'classification/rgc', 'onlyL', ...
    'classifySVM_rgcResponse_Cones2RGC5_absorptionrate.mat'));

for ratio = 1:5
    % Load WITH stim template classifier results
    wTemplate=load(fullfile(baseFolder, 'data', expName, 'classification/rgc',  'stimTemplate',subFolder, ...
        sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_defaultnophaseshift_%s.mat',ratio,eccen,subFolder)));
    
    performance.templateEnergy = wTemplate.P_svmEnergy;
    performance.templateLinear = wTemplate.P_svmLinear;
    performance.noTemplate = woTemplate.P_svm(ratio,:);
    
    logzero = 5.^-6;
    
    subplot(3,1,1);
    plot(wTemplate.expParams.contrastLevels, performance.templateEnergy, 'o-', 'color',cmap(ratio,:),  'lineWidth',2); hold on;
    
    subplot(3,1,2);
    plot(wTemplate.expParams.contrastLevels, performance.templateLinear, 'o-', 'color',cmap(ratio,:),'lineWidth',2); hold on;
    
    
    subplot(3,1,3);
    plot(woTemplate.expParams.contrastLevels, performance.noTemplate, 'o-', 'color',cmap(ratio,:),'lineWidth',2); hold on;
    
end

% Add labels, titles and set axis limits
for r = 1:5; labels{r} = sprintf('%1.2f RGC : 1 cone',2/r); end

subplot(311);
legend(labels, 'Location', 'NorthEast')
legend boxoff
for ratio = 1:5
    plot(logzero,  performance.templateEnergy(1), 'o', 'color',cmap(ratio,:));
end
ylim([40 100]); xlim(10.^[-4.5 -1])
set(gca, 'XScale', 'log')
xlabel('Stimulus contrast (fraction)')
ylabel('Accuracy (%)')
set(gca, 'TickDir', 'out', 'FontSize',15)
title({sprintf('SVM-Energy Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})

subplot(312);
legend(labels, 'Location', 'NorthEast')
legend boxoff
for ratio = 1:5
    plot(logzero,  performance.templateLinear(1), 'o', 'color',cmap(ratio,:));
end
ylim([40 100]); xlim(10.^[-4.5 -1]);
set(gca, 'XScale', 'log')
xlabel('Stimulus contrast (fraction)')
ylabel('Accuracy (%)')
set(gca, 'TickDir', 'out', 'FontSize',15)
title({sprintf('SVM-Linear Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})

subplot(313)
legend(labels, 'Location', 'NorthEast')
legend boxoff
for ratio = 1:5
    plot(logzero,  performance.noTemplate(1), 'o', 'color',cmap(ratio,:));
end
ylim([40 100]); xlim(10.^[-4.5 -1])
set(gca, 'XScale', 'log')
xlabel('Stimulus contrast (fraction)')
ylabel('Accuracy (%)')
set(gca, 'TickDir', 'out', 'FontSize',15)
title({sprintf('SVM-Fourier Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})



%% Try for cone density RGC simulations:
% Eccen: 0-40 deg, LMS cone mosaics, drift and microsaccades, 2 phases
% within an oriented Gabor class, and Poisson noise at cone level.

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'conedensity';
seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
ratio      = 1; % cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell)
eccen      = 5; % equal to 4.5 deg eccen

[PEnergy, PLinear] = linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen);

%% Visualize differences in psychometric functions
cmap = jet(13);
figure(2); clf; set(gcf, 'Position',[ 332    39   618   759], 'color','w'); hold all;
for eccen = 1:13
    for ratio = 1:5
        % Load WITH stim template classifier results
        wTemplate=load(fullfile(baseFolder, 'data', expName, 'classification','rgc', 'meanPoissonPadded', 'stimTemplate',subFolder, ...
            sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_%s.mat',ratio,eccen,subFolder)));
        
        % Load WITHOUT stim template classifier results
        woTemplate=load(fullfile(baseFolder, 'data', expName, 'classification/rgc', 'withPaddingBeforeConvolution', ...
            [subFolder '_meanPoissonPadded'], ...
            sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_%s_meanPoissonPadded.mat',ratio, eccen,subFolder)));
        
        performance.templateEnergy = wTemplate.P_svmEnergy;
        performance.templateLinear = wTemplate.P_svmLinear;
        performance.noTemplate = woTemplate.P_svm;
        
        logzero = 1.^-4;
        
        if (eccen >= 10) && (ratio ==5)
            wTemplate.expParams.contrastLevels = wTemplate.expParams.contrastLevelsPC;
            woTemplate.expParams.contrastLevels = woTemplate.expParams.contrastLevelsPC;
        end
        
        subplot(3,5,ratio); 
        plot(wTemplate.expParams.contrastLevels, performance.templateEnergy, 'o-', 'color',cmap(eccen,:),  'lineWidth',2); hold all;
        plot(logzero,  performance.templateEnergy(1), 'o', 'color',cmap(eccen,:));
        ylim([40 100]); xlim(10.^[-3 -1])
        set(gca, 'XScale', 'log')
        title('SVM-Energy')
    
        subplot(3,5,ratio+5);
        plot(wTemplate.expParams.contrastLevels, performance.templateLinear, 'o-', 'color',cmap(eccen,:),'lineWidth',2); hold all;
        plot(logzero,  performance.templateLinear(1), 'o', 'color',cmap(eccen,:));
         ylim([40 100]); xlim(10.^[-3 -1])
        set(gca, 'XScale', 'log')
        title('SVM-Linear')
    
        subplot(3,5,ratio+10); 
        plot(woTemplate.expParams.contrastLevels, performance.noTemplate, 'o-', 'color',cmap(eccen,:),'lineWidth',2); hold all;
        plot(logzero,  performance.noTemplate(1), 'o', 'color',cmap(eccen,:));
        ylim([40 100]); xlim(10.^[-3 -1])
        set(gca, 'XScale', 'log');
        title('SVM-Fourier')
    
    end
end

%% Add labels, titles and set axis limits
for r = 1:5; labels{r} = sprintf('%1.2f RGC : 1 cone',2/r); end

for ratio = [5,10,15]
    subplot(2,5,ratio);
    legend(labels, 'Location', 'BestOutside')
    legend boxoff
end


   
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    title({sprintf('SVM-Energy Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})

   

    ylim([40 100]); xlim(10.^[-4.5 -1]);
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    title({sprintf('SVM-Linear Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})

    ylim([40 100]); xlim(10.^[-4.5 -1])
    set(gca, 'XScale', 'log')
    xlabel('Stimulus contrast (fraction)')
    ylabel('Accuracy (%)')
    set(gca, 'TickDir', 'out', 'FontSize',15)
    title({sprintf('SVM-Fourier Classifier performance for RGC responses'),'4.5 deg eccen, L-only, no eye movements, no stimulus phase shifts'})

