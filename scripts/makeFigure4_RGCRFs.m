function [] = makeFigure4_RGCRFs()

baseFolder  = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName     = 'default';
ratios      = 1:5; % corresponding to cone:RGC ratio's 2:1, 1:1, 0.67:1, 0.5:1, 0.4:1
labelsRatio = sprintfc('cone2RGC = %1.1f:1', 2./ratios);
colors      = parula(6);


%% Visualize visual field of 1D DoG filter

fH1 = figure(1); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH2 = figure(2); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;

for r = ratios

    % Load data
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
        
    % Get size, midpoint and quarterpoints of array
    sz = size(DoGfilter);
    midpoint = ceil(sz(1)/2);
    quarterpoint = midpoint/2;
    
    % Build cone and corresponding rgc array
    [X,Y] = meshgrid(1:rgcParams.cRows,1:rgcParams.cCols);
    conearray = zeros(rgcParams.cCols, rgcParams.cRows);
    rowIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cRows;
    colIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cCols;
    rgcarray = conearray;
    rgcarray(rowIndices, colIndices) = 1;
    
    % X-axis for 1D DoG visual field
    x = linspace(1,sz(1),numel(DoGfilter(midpoint,:)));
    
    % 1D DOG visual field
    figure(1)
    subplot(1,length(ratios),r);
    plot([x(1) x(end)], [0 0], 'k'); hold on;
    plot(x,DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    plot(x+r, DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    plot(x+r, DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    plot(x-r, DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    plot(x-r, DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)

    set(gca, 'YLim', [-0.05,1],'XLim', (midpoint + [-10, 10]), 'TickDir', 'out', 'FontSize', 15, 'XTick', [midpoint-10,midpoint,midpoint+10], ...
        'XTickLabel', {'-0.25', '0', '0.25'});
    xlabel('x position (deg)')
    ylabel('modulation (a.u.)')
    title(labelsRatio(r))
    axis square; box off;
    
    
    % 2D Array visual field
    figure(2)
    subplot(1,length(ratios),r);
    plot(X,Y, '.','color', [0.5,0.5,0.5], 'MarkerSize',9); hold on;
    plot(X(logical(rgcarray)),Y(logical(rgcarray)), 'r.');
    axis square;box off; 
    title(labelsRatio{r});
    xlabel('x position (deg)');
    ylabel('y position (deg)');
    set(gca, 'XLim', (midpoint + [-10, 10]),'YLim', (midpoint + [-10, 10]), ...
        'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [midpoint-10,midpoint,midpoint+10], ...
        'XTickLabel', {'-0.25', '0', '0.25'}, ...
        'YTick', [midpoint-10,midpoint,midpoint+10], ...
        'YTickLabel', {'-0.25', '0', '0.25'});
    
end
    
hgexport(fH1, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs'))
print(fH1, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs'), '-dpng')
hgexport(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'))
print(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'), '-dpng')


%% Visualize FFT of 1D DoG filter 

fH3 = figure(3); set(gcf, 'Position', [786,135,560,420], 'color', 'w'); clf; hold all;

for r = ratios
    
    % Load data
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
    
    G_2D = abs(fft2(DoGfilter./sum(abs(DoGfilter(:)))));
    
    % Visualize  1D representation
    ncones = size(G_2D,1);
    fs     = (0:ncones-1)/2;
    G      = fftshift(G_2D);
    
    midpoint = ceil(size(G,1)/2);
    G_1D  = G(:,midpoint);

    plot(fs-(max(fs)/2), G_1D, 'color', colors(r,:), 'LineWidth', 4); hold on;
end

    plot(rgcParams.stimSF * [1 1], [0 .8], 'r:', 'LineWidth', 4);
    xlabel('Spatial frequency (cycles/deg)'); 
    ylabel('Normalized Modulation (a.u.)'); 
    legend(labelsRatio, 'Location', 'Best'); legend boxoff;
    set(gca, 'YLim', [0 .8],'XLim', [0 max(fs)/2], ...
            'TickDir', 'out', 'FontSize', 15);
    axis square; box off;

hgexport(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_DoGFFT'))
print(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_DoGFFT'), '-dpng')
