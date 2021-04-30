function makeFigure6_3DThresholdDensityRatio()
% Function to make Figure 6 of the manuscript:
%   Radial asymmetries around the visual field: From retina to cortex to 
%   behavior. By Kupers, Benson, Carrasco, Winawer.
%    JOURNAL. DOI.

% This function requires you to have psychometric functions of the 
% simulated data. You can run this with:
% 
% for ratio = 1:5
%     plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',ratio, 'plotAvg',true, 'meanPoissonPaddingFlag', true);
% end

%% 0. Define params and folders
meanPoissonPaddingFlag = true;
stimTemplateFlag       = false;
plotRGC2Cones          = true; % if false, plot cone2rgcs
saveFigs               = true;

nrSimConeDensities     = 13;
lowessSpan             = 0.4; % how smooth should 3D fit be
c2rgc                  = (1:5).^2; % Cone 2 RGC ratios
expName                = 'conedensity';
colors                 = parula(length(c2rgc)+1);

if plotRGC2Cones
    ratiosToPlot       = 2./c2rgc; % RGC 2 Cone ratios
    labels             = sprintfc('RGC:cone = %1.2f:1.0', ratiosToPlot);
else
    ratiosToPlot       = c2rgc;
    labels             = sprintfc('Cone:RGC = 1.0:%1.2f', ratiosToPlot);
end


% Change folder names if using mean Poisson padded cone data
if meanPoissonPaddingFlag
    extraSubFolder = 'meanPoissonPadded';
else
    extraSubFolder = 'noPaddingBeforeConvolution';
end

if stimTemplateFlag
    subFolder = 'SVM-Energy';
else
    subFolder = 'SVM-Fourier';
end

% Folders
baseFolder   = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model';
dataFolder   = fullfile(baseFolder,'data',expName,'thresholds','rgc',extraSubFolder,subFolder,'average');
if plotRGC2Cones
    figureFolder = fullfile(baseFolder,'figures','surface3D_fixedRatios_rgc2c', expName, extraSubFolder, subFolder,'average');
else
    figureFolder = fullfile(baseFolder,'figures','surface3D_fixedRatios_c2rgc', expName, extraSubFolder, subFolder,'average');
end
if (saveFigs) && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

%% 1. Load thresholds for each cone2mrgc ratio and fit function to data

% preallocate space
allData   = NaN(length(c2rgc), nrSimConeDensities);
errThresh = allData;

% Loop over ratios
for r = 1:length(c2rgc)
    
    % Load simulated contrast thresholds
    load(fullfile(dataFolder, sprintf('cThresholds_ratio%d_average.mat', r)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');
    allData(r,:) = [reshape(cell2mat(fit.ctrthresh),[],length(fit.ctrthresh))].*100;
    
    % Load variance in simulated contrast thresholds computed by
    % bootstrapping simulation iterations starting with different rng seeds
    load(fullfile(dataFolder, sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptions_13_conedensity', r)), 'varThresh');
    errThresh(r,:) = varThresh.*100;
    
%     contrastThresh = cell2mat(fit.ctrthresh).*100; 
    coneDensities  = xThresh';

end
% clean up 
clear fit xThresh fitToPlot dataToPlot;

%% Construct 3D fit

% Find NaNs
idx = isfinite(allData');

% Make 2D grid
[x,y] = meshgrid(log10(c2rgc),log10(coneDensities));
z = log10(allData');

% Fit it!
[meshFit, gof] = fit([x(idx) y(idx)], z(idx), 'lowess','span',lowessSpan);
R2_mesh = gof.rsquare;

% Extract single lines for separate ratio's
singleRatioMeshfit = meshFit(x,y);
% R2 = gof.rsquare;

%% ---------------------------------------------------------------
% ---------------- Visualize single ratio results -------------- %
% ----------------------------------------------------------------

% Plot threshold fits vs density data for each ratio
fH1 = figure(1); set(gcf, 'Color', 'w', 'Position', [79 39 1522 759]); clf; hold all;
if stimTemplateFlag
    ylimit = [0.5 100]; 
    yticks = [1 10 100];
else
    ylimit = [0.5, 20];
    yticks = [1 10];
end
for ii = 1:length(c2rgc)
    subplot(2,3,ii); hold on;
    errorbar(coneDensities, allData(ii,:), errThresh(ii,:), 'Color', 'k', 'LineStyle','none', 'LineWidth', 0.5, 'CapSize',2);
    scatter(coneDensities,  allData(ii,:), 30, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',0.5)
    plot(coneDensities, 10.^singleRatioMeshfit(:,ii)', 'r-', 'LineWidth',2)
    
    set(gca, 'XScale', 'log','YScale','log', 'FontSize', 14, 'LineWidth', 0.5, 'TickDir', 'out', ...
        'YTick',[1 10], 'YTickLabel',yticks)
    ylim(ylimit);
    title(labels{ii});
    xlabel('Cone array density (cells/deg^2)');
    ylabel('Contrast thresholds (%)');
    box off; grid on;
end

% Plot lines in additional subplot
subplot(2,3,6); hold all;
for ii = 1:length(c2rgc)
    plot(coneDensities, 10.^singleRatioMeshfit(:,ii)', '-', 'LineWidth',4, 'Color', colors(ii,:))
end

h = findobj(gca,'Type','line');
legend([h(end:-1:1)],labels, 'Location','Best');
legend boxoff;
set(gca, 'XScale', 'log','YScale','log', 'FontSize', 14, 'LineWidth', 0.5, 'TickDir', 'out', ...
    'YTick',ylimit, 'YTickLabel',yticks)
ylim([0.5, 20]); xlim([10.^2,10.^5]);
title('All 5 mRGC : cones ratios')
xlabel('Cone array density (cells/deg^2)');
ylabel('Contrast thresholds (%)');
box off; grid on;

if saveFigs
    hgexport(fH1,fullfile(figureFolder, sprintf('ContrastThreshold-vs-Density_Lowess%1.2f.eps', lowessSpan)))
    savefig(fH1, fullfile(figureFolder, sprintf('ContrastThreshold-vs-Density_Lowess%1.2f.fig', lowessSpan)))
    saveas(fH1, fullfile(figureFolder, sprintf('ContrastThreshold-vs-Density_Lowess%1.2f.png', lowessSpan)))
end


%% ---------------------------------------------------------------------
% ---------------- Visualize mesh of all ratios results -------------- %
% ----------------------------------------------------------------------

fH3 = figure(3); clf; set(gcf, 'Position', [782 44 881 756], 'Color', 'w'); clf;
ax = plot(meshFit,[x(idx) y(idx)], z(idx));
ax(1).FaceLighting = 'gouraud';
ax(1).FaceColor = [1 1 1];
ax(2).Marker = 'none';

% Add labels
if plotRGC2Cones
    xlabel('mRGC:cone ratio', 'FontSize', 20)
else
    xlabel('Cone:mRGC ratio', 'FontSize', 20)
end
ylabel('Cone density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Interpolated fit to data: Effect of RGC filtering on contrast threshold')

% Make plot pretty
set(gca, 'ZLim',[log10(0.7) 1.1],'FontSize', 20, 'LineWidth', 2, ...
    'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear', ...
    'XTick', log10(c2rgc), 'XTickLabel', sprintfc('%1.2f', ratiosToPlot), 'XDir','reverse',...'YDir','reverse',...
    'YLim',[2 5],'YTick',[2:1:5],'YTickLabel',{'10^2','10^3','10^4','10^5'},...
    'ZTick', log10([0.7:0.1:1,2:1:10]), 'ZTickLabel',[0.7:0.1:1,2:1:10],...
    'TickDir', 'out','View',[-134.4000   11.2000]);
grid on;
set(gca, 'GridAlpha', .1, 'ZMinorGrid', 'off', 'YMinorGrid', 'off', 'XMinorGrid', 'off')
set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1])
axis square; material shiny;

% Save figure if  requested
if saveFigs
    hgexport(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view1.eps',lowessSpan)))
    savefig(fH3, fullfile(figureFolder,  sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view1.fig',lowessSpan)))
    print(fH3, fullfile(figureFolder,  sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view1.png',lowessSpan)), '-dpng')
end

% Flip view and save again
% set(gca, 'View', [-134.4000   11.2000])

% if saveFigs
%     hgexport(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view2.eps',lowessSpan)))
%     savefig(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view2.fig',lowessSpan)))
%     print(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_view2.png',lowessSpan)), '-dpng')
% end

%%  Predict contrast thresholds for given mRGC density

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

if plotRGC2Cones
    % Compute cone:RGC ratio
    rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;
    ratioAtIdx   = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index
    
    % Get cone density at chosen eccentricity for each meridian
    observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;
    
    % Check: should be equal to curcio data
    isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));
    
    % take reciprocal for plotting -- meshfit expects cone2rgc ratio
    ratioAtIdx = (1./ratioAtIdx);
else
    % Compute RGC:cone ratio
    cone2RGCRatio = coneDensityDeg2PerMeridian./watson2015.mRGCRFDensityPerDeg2;
    ratioAtIdx = cone2RGCRatio(:,idxEccen);
    
    % Get cone density at chosen eccentricity for each meridian
    observedConesAtEccen = ratioAtIdx.*watson2015.mRGCRFDensityPerDeg2(:,idxEccen);
    
    % Check: should be equal to curcio data
    isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));
end

% Find contrast threshold data for all meridians: Nasal, Superior,temporal, inferior
predictedContrastThreshold = meshFit(log10(ratioAtIdx),log10(observedConesAtEccen));

%% Plot observed/biological variations in cone:mRGC ratios on mesh
figure(fH3); hold all;
colorsRetina = {'r', 'b', 'g', 'k'};
zlift        = [0.01, 0.01, 0.01, 0.01]; % lift markers a tiny bit for visibility

for jj = 1:4
    scatter3(log10(ratioAtIdx(jj)),log10(observedConesAtEccen(jj)),predictedContrastThreshold(jj)+zlift(jj), 300, 'MarkerFaceColor', colorsRetina{jj}, 'MarkerEdgeColor','k', 'LineWidth',0.1, 'Marker', 'p')
end

% Save figure 3 again with dots
if saveFigs
    hgexport(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view2.eps', lowessSpan)))
    savefig(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view2.fig', lowessSpan)))
    print(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view2.png', lowessSpan)), '-dpng')
end

% Ratio left / cone density right
set(gca, 'View', [-45.6000   11.2000])
set(gca, 'ydir', 'reverse')

if saveFigs
    hgexport(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view1.eps',lowessSpan)))
    savefig(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view1.fig',lowessSpan)))
    print(fH3, fullfile(figureFolder, sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view1.png',lowessSpan)), '-dpng')
end

%% Make gif of rotating mesh
%
% filename = fullfile(figureFolder,'fullModelMesh.gif');
%
% x = -160:5:-95;
% angles = [x, fliplr(x)];
%
% axis manual; axis vis3d;
% for n = 1:length(angles)
%     set(gca, 'View', [angles(n), 10.4000]);
%
%     drawnow; pause(0.1);
%
%     % Capture the plot as an image
%     frame = getframe(fH3);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if n == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
%     end
% end


