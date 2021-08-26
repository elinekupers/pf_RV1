function fH = makeFigure_3DSurfaceMesh_LateNoiseRGCModel(...
            rgc3D, ratioAtIdx, observedConesAtEccen, ...
            predictedContrastThreshold, expName, figurePth, saveFig)

% Get downsample factors as x-axis
downsampleFactors = 2./(1:5).^2; % RGC:cone downsample ratios for 2D arrays
xticks = fliplr(downsampleFactors); % get x axis range and xtick labels
for ii = 1:length(xticks)
    xlabels_downsample{ii} = sprintf('%1.2f', xticks(ii));
end

fH = figure(12); clf; set(gcf, 'Position', [782 44 881 756], 'Color', 'w', ...
    'NumberTitle', 'off', 'Name', sprintf('3D surface: Contrast threshold vs downsample factor vs cone density: %s', expName));
ax = plot(rgc3D.meshFit,[rgc3D.X(:) rgc3D.Y(:)], rgc3D.Z(:));
ax(1).FaceLighting = 'gouraud';
ax(1).FaceColor = [1 1 1];
ax(2).Marker = 'none';

% Add labels
xlabel('mRGC : cone ratio', 'FontSize', 20)
ylabel('Cone density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Effect of late noise RGC model on contrast threshold')

zticks = [0.001:0.001:0.01, ...
    0.02:0.01:0.1];
ztick_labels = cell(size(zticks));
printTick = [0.01, 0.1];
for ii = 1:length(printTick)
    ztick_labels{zticks==printTick(ii)} = sprintf('%2.0f',printTick(ii)*100);
end

% Make plot pretty
set(gca, 'ZLim',[-2.5 -1],'FontSize', 20, 'LineWidth', 2, ...
    'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear', ...
    'XTick', log10(fliplr(downsampleFactors)), ...
    'XTickLabel', xlabels_downsample, 'XDir','reverse',...'YDir','reverse',...
    'YLim',[2 4],'YTick',[2:1:4],'YTickLabel',{'10^2','10^3','10^4',},...
    'ZTick', log10(zticks), 'ZTickLabel',ztick_labels,...
    'TickDir', 'out','View',[-134.4000   11.2000]);
grid on;
set(gca, 'GridAlpha', .2, 'ZMinorGrid', 'off', 'YMinorGrid', 'off', 'XMinorGrid', 'off')
set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1])
axis square; material shiny;

% Plot observed biological variations in cone:mRGC ratios on mesh
hold all;
colorsRetina = {'r', 'b', 'g', 'k'}; % nasal, superior, temporal, inferior retina
zlift        = 2*[0.01, 0.01, 0.01, 0.01]; % lift markers a tiny bit for visibility

for jj = 1:4
    scatter3(log10(ratioAtIdx(jj)),log10(observedConesAtEccen(jj)), ...
        log10(predictedContrastThreshold(jj))+zlift(jj), 300, ...
        'MarkerFaceColor', colorsRetina{jj}, 'MarkerEdgeColor','k', ...
        'LineWidth',0.1, 'Marker', 'p')
end

% Save figure
if saveFig
    fName = sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFit%s_withDots_view1',rgc3D.whichfit);
    hgexport(fH, fullfile(figurePth, [fName '.eps']))
    savefig(fH, fullfile(figurePth, [fName '.fig']))
    print(fH, fullfile(figurePth, [fName '.png']), '-dpng')
end

% Rotate view: Ratio left / cone density right
set(gca, 'View', [-45.6000   11.2000])
set(gca, 'ydir', 'reverse')

% Save figure again
if saveFig
    fName = sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFit%s_withDots_view2',rgc3D.whichfit);
    hgexport(fH, fullfile(figurePth, [fName '.eps']))
    savefig(fH, fullfile(figurePth, [fName '.fig']))
    print(fH, fullfile(figurePth, [fName '.png']), '-dpng')
end