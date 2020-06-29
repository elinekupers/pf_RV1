% makeFigure4_ModelOverview

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName = 'default';
ratios = 1:5;
labelsRatio = sprintfc('cone2RGC = %1.1f:1', 2./ratios);
colors = parula(6);

fH1 = figure(1); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH2 = figure(2); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH11 = figure(11); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;

%% Visualize filter and subsample grid
for r = ratios

    % Load data
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
    
    labelsFilter = sprintf('ctr:sur = %1.2f:%1.2f', rgcParams.DoG.kc,rgcParams.DoG.ks);
    
    % Get array
    sz = size(DoGfilter);
    midpoint = ceil(sz(1)/2);
    quarterpoint = midpoint/2;
    
    [X,Y] = meshgrid(1:rgcParams.cRows,1:rgcParams.cCols);
    conearray = zeros(rgcParams.cCols, rgcParams.cRows);
    rowIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cRows;
    colIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cCols;
    rgcarray = conearray;
    rgcarray(rowIndices, colIndices) = 1;
    
    %DOG visual field
    figure(1)
    subplot(1,length(ratios),r);
    plot([0,sz(1)], [0 0], 'k'); hold on;
    plot(DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    xlabel('x position (deg)')
    ylabel('modulation (a.u.)')
    title(labelsFilter)
    set(gca, 'YLim', [-0.05,1],'XLim', (midpoint + [-10, 10]), 'TickDir', 'out', 'FontSize', 15, 'XTick', [midpoint-10,midpoint,midpoint+10], ...
        'XTickLabel', {'-0.25', '0', '0.25'});
    axis square; box off;
    
    figure(11)
    subplot(1,length(ratios),r);
    imagesc(DoGfilter); colormap gray; colorbar;
    xlabel('x position (deg)')
    ylabel('x position (deg)')
    title(labelsFilter)
    set(gca, 'CLim', [-0.05, 1], 'YLim', [0 sz(1)], 'XLim', [0 sz(1)], 'TickDir', 'out', 'FontSize', 15, 'XTick', [midpoint-10,midpoint,midpoint+10], ...
        'XTickLabel', {'-0.25', '0', '0.25'});
    axis square; box off;
    
    %Array visual field
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
hgexport(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'))
hgexport(fH11, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs_2D_VisualField'))

print(fH1, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs'), '-dpng')
print(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'), '-dpng')
print(fH11, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs_2D_VisualField'))


%% All DoG filters in one figure
fH3 = figure(3); set(gcf, 'Position', [786,135,560,420], 'color', 'w'); clf; hold all;

xlabel('x position (deg)')
ylabel('modulation (a.u.)')
title('DoG filters')
set(gca,'YLim', [-0.05,1], 'XLim', [0, sz(1)], 'TickDir', 'out', 'FontSize', 15, 'XTick', [0, midpoint, sz(1)], ...
    'XTickLabel', {'-1', '0', '1'});
axis square; box off;

for r = ratios
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
    plot(DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2); hold all;
end
legend(labelsRatio, 'Location','Best'); legend boxoff;

plot([0,sz(1)], [0 0], 'k'); hold on;


hgexport(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs_allInOne'))
print(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs_allInOne'), '-dpng')

%% 2D RGC RESPONSE -- FILTERED ONLY 

fH4 = figure(4); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
for r = ratios
    
    fNameFiltResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', 'filteredOnly', sprintf('filteredConeCurrent_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltResponse), 'filteredConeCurrent');
    midpoint = ceil(size(filteredConeCurrent,1)/2);
    
    subplot(1,length(ratios),r); cla
    imagesc(filteredConeCurrent); colormap gray; hold on;
    plot([0, size(filteredConeCurrent,1)], [midpoint, midpoint], 'r:', 'LineWidth',4)
    colorbar;
    xlabel('x position (deg)')
    ylabel('y position (deg)')
    title(labelsRatio{r});
    set(gca, ...
        'CLim', [0 max(filteredConeCurrent(:))])
    set(gca, 'YLim', [0 sz(1)],'XLim', [0 sz(1)], ...
        'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [0, midpoint, size(filteredConeCurrent,1)], ...
        'XTickLabel', {'-1', '0', '1'}, ...
        'YTick', [0, midpoint, size(filteredConeCurrent,2)], ...
        'YTickLabel', {'-1', '0', '1'});
    axis square; box off;
end

hgexport(fH4, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseFilteredOnly'))
print(fH4, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseFilteredOnly'), '-dpng')


%% 2D RGC RESPONSE FILTERED + SUBSAMPLE 
fH5 = figure(5); set(gcf, 'Position', [33,42,1639,763], 'color', 'w'); clf; hold all;

for r = ratios
    
    fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', sprintf('ratio%d',r), sprintf('rgcResponse_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltSubResponse));

    % Use selected time points only
    filteredSubsampled = squeeze(mean(rgcResponse(:,:,:,[1:3],1),1))./2;
    midpoint = ceil(size(filteredSubsampled,1)/2);

    for ii = 1:3
        subplot(3,length(ratios),((ii-1)*length(ratios))+r);
        imagesc(squeeze(filteredSubsampled(:,:,ii))); colormap gray; hold on;
        
        plot([0, size(filteredConeCurrent,1)], [midpoint, midpoint], 'r:', 'LineWidth',4)
        colorbar;
        xlabel('x position (deg)')
        ylabel('y position (deg)')
        title(labelsRatio{r});
        set(gca, ...
            'CLim', [0 max(filteredSubsampled(:))], 'TickDir', 'out', 'FontSize', 15)
        axis square; box off;
    end
end  
    
hgexport(fH5, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseSubsampled'))
print(fH5, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseSubsampled'), '-dpng')


%% 1D RGC RESPONSE FILTERED + SUBSAMPLE 
fH6 = figure(6); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;

for r = ratios
    
    fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', sprintf('ratio%d',r), sprintf('rgcResponse_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltSubResponse));

    % Use selected time points only
    midpoint = ceil(size(rgcResponse,2)/2);

    filteredSubsampled_1D = squeeze(mean(rgcResponse(:,:,midpoint,1,1),1))./2;
    
    
    subplot(1,length(ratios),r);  hold on;
    plot(filteredSubsampled_1D, 'k', 'LineWidth', 2); 
    xlabel('x position (deg)', 'FontSize', 15)
    ylabel('RGC response (a.u.)', 'FontSize', 15)
    title(labelsRatio{r});
    set(gca, 'YLim', [0 max(filteredSubsampled_1D)],'XLim', [0 length(filteredSubsampled_1D)], ...
            'TickDir', 'out', 'FontSize', 15, ...
            'XTick', [0, midpoint, length(filteredSubsampled_1D)], ...
            'XTickLabel', {'-x', '0', 'x'});
    axis square; box off;
    
end 

hgexport(fH6, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_rgcResponseSubsampled'))
print(fH6, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_rgcResponseSubsampled'), '-dpng')


%% 1D DOG FFT

fH7 = figure(7); set(gcf, 'Position', [786,135,560,420], 'color', 'w'); clf; hold all;

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

hgexport(fH7, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_DoGFFT'))
print(fH7, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_DoGFFT'), '-dpng')

    
%% 2D DOG FFT
    
fH8 = figure(8); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
for r = ratios
    
    % Load data
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));

    G_2D = abs(fft2(DoGfilter./sum(abs(DoGfilter(:)))));
    midpoint = ceil(size(G_2D,1)/2);
    G_1D  = G_2D(:,midpoint);
   
    % get stim sf and get coords to draw circle
    sfOfStim = rgcParams.stimSF * rgcParams.fov;

    th = 0:pi/50:2*pi;
    xvals = sfOfStim * cos(th) + midpoint;
    yvals = sfOfStim * sin(th) + midpoint;
    
    subplot(1,length(ratios),r); cla; hold on;
    imagesc(fftshift(G_2D)); colormap gray;
    
    % draw center cross
    plot([midpoint, midpoint], [1,size(G_2D,1)], ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);
    plot([1,size(G_2D,1)],[midpoint, midpoint],  ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);

    % draw stim sf cross
    plot(xvals, yvals, 'r:', 'lineWidth',3);
    colorbar; 
    set(gca, 'CLim', [0 0.8], 'XLim', [0 size(G_2D,1)])
    xlabel('cycles/deg');
    ylabel('cycles/deg');
    title(labelsRatio{r})
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
       'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'XTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)}, ...
        'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'YTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)});
    axis square; box off;
    
end

hgexport(fH8, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_DoGFFT'))
print(fH8, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_DoGFFT'), '-dpng')


%% 1D and 2D RGC RESPONSE FFT

fH9 = figure(9); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH10 = figure(10); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;

for r = ratios
    
    fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', sprintf('ratio%d',r), sprintf('rgcResponse_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltSubResponse));

    % Use selected time points only
    filteredSubsampled = squeeze(mean(rgcResponse(:,:,:,1,1),1))./2;
    filteredSubsampled_amps_2D  = abs(fft2(filteredSubsampled));
    filteredSubsampled_amps_2D(1,1) = NaN;
    filteredSubsampled_amps_2D = fftshift(filteredSubsampled_amps_2D);
%     filteredSubsampled_amps_2D_norm = fftshift(filteredSubsampled_amps_2D)./max(filteredSubsampled_amps_2D(:));
    
    midpoint = 1+ceil(size(filteredSubsampled_amps_2D,1)/2);  
    quarterpoint_rsp = midpoint/2;
    
    figure(9);
    subplot(1,length(ratios),r); cla
    imagesc(filteredSubsampled_amps_2D); colormap gray;
    colorbar; 
    set(gca, ...
         'CLim', [0 max(filteredSubsampled_amps_2D(:))])
    xlabel('sf (cycles/deg)')
    ylabel('sf (cycles/deg)')
    title(labelsRatio{r})
      
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'XTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)}, ...
        'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'YTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)});
    axis square; box off;



    
    filteredSubsampled_amps_1D  = filteredSubsampled_amps_2D(:,midpoint);
    ncones = length(filteredSubsampled_amps_1D);
    x      = (0:ncones-1)/ncones*rgcParams.fov; % degrees
    fs     = (0:ncones-1)/2;
    
    figure(10);
    subplot(1,length(ratios),r); cla; hold on;
    plot(fs-max(fs)/2, filteredSubsampled_amps_1D, 'k', 'LineWidth', 2); 
    xlim([1 max(fs)/2])
    plot(rgcParams.stimSF * [1 1], [0 1], 'r:', 'LineWidth', 4);
    xlabel('sf (cycles/deg)')
    ylabel('modulation (a.u.)')
    title(labelsRatio{r})
    set(gca, ...
        'TickDir', 'out', 'FontSize', 15, 'XTick', [1, 4, max(fs)/2], ...
        'XTickLabel', {'1', '4', num2str(max(fs)/2)});
    axis square; box off;
%     
 
end

hgexport(fH9, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseFFT'))
print(fH9, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '2D_rgcResponseFFT'), '-dpng')
hgexport(fH10, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_rgcResponseFFT'))
print(fH10, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', '1D_rgcResponseFFT'), '-dpng')



