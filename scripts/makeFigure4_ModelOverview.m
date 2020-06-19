% makeFigure4_ModelOverview

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName = 'default';
ratios = 1:5;
labelsRatio = sprintfc('cone2RGC = %1.1f:1.0', 2./ratios);
colors = parula(6);

fH1 = figure(1); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH2 = figure(2); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH3 = figure(3); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH4 = figure(4); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH5 = figure(5); set(gcf, 'Position', [786,135,560,420], 'color', 'w'); clf; hold all;
fH6 = figure(6); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH7 = figure(7); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH8 = figure(8); set(gcf, 'Position', [33,719,1639,234], 'color', 'w'); clf; hold all;
fH9 = figure(9); set(gcf, 'Position', [786,135,560,420], 'color', 'w'); clf; hold all;

for r = ratios

    % Load data
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
    
    labelsFilter = sprintf('ctr:sur = %1.2f:%1.2f', rgcParams.DoG.kc,rgcParams.DoG.ks);
    
    fNameFiltResponse = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('filteredConeCurrent_Cones2RGC%d_contrast0.1000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltResponse));
   
    fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d_contrast0.1000_eccen4.50_absorptionrate.mat',r));
    load(fullfile(fNameFiltSubResponse));
    
    % Use selected time points only
    timepointsToAverage = 1:28;
    filteredSubsampled = squeeze(mean(mean(rgcResponse(:,:,:, timepointsToAverage,1),1),4))./2;
    
    % Get array
    sz = size(DoGfilter);
    midpoint = ceil(sz(1)/2);
    [X,Y] = meshgrid(1:rgcParams.cRows,1:rgcParams.cCols);
    conearray = zeros(rgcParams.cCols, rgcParams.cRows);
    rowIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cRows;
    colIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cCols;
    rgcarray = conearray;
    rgcarray(rowIndices, colIndices) = 1;
    
    % Stim and noise in FFT
    filteredOnly_amps_1D  = abs(fft(filteredConeCurrent(:,midpoint)));
    
    filteredSubsampled_amps_2D  = abs(fft2(filteredSubsampled));
    midpoint_rsp = ceil(size(filteredSubsampled_amps_2D,1)/2);
    filteredSubsampled_amps_1D = filteredSubsampled_amps_2D(:,midpoint_rsp);
    
    % Get FFT amplitudes of Gaussian filter 
    G_2D = abs(fft2(DoGfilter./sum(abs(DoGfilter(:)))));
    G_1D  = G_2D(:,midpoint);

    quarterpoint = midpoint/2;
    quarterpoint_rsp = midpoint_rsp/2;

    % get stim sf and get coords to draw circle
    sfOfStim = rgcParams.stimSF * rgcParams.fov;

    th = 0:pi/50:2*pi;
    xvals = sfOfStim * cos(th) + midpoint;
    yvals = sfOfStim * sin(th) + midpoint;
    
    %% DOG visual field
   
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
    
    
    figure(9); hold all;
    plot(DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2)
    
    
    %% 1D RGC RESPONSE
    figure(3)
    subplot(1,length(ratios),r); cla
    imagesc(filteredSubsampled); colormap gray; hold on;
    plot([0, size(filteredSubsampled,1)], [midpoint_rsp, midpoint_rsp], 'r:', 'LineWidth',4)
    colorbar;
    xlabel('x position (deg)')
    ylabel('y position (deg)')
    title(labelsFilter);
    set(gca, ...
        'CLim', [0 max(filteredSubsampled(:))])
%     set(gca, 'YLim', [0 sz(1)],'XLim', [0 sz(1)], ...
%         'TickDir', 'out', 'FontSize', 15, ...
%         'XTick', [0, midpoint_rsp, size(filteredSubsampled,1)], ...
%         'XTickLabel', {'-1', '0', '1'}, ...
%         'YTick', [0, midpoint_rsp, size(filteredSubsampled,2)], ...
%         'YTickLabel', {'-1', '0', '1'});
    axis square; box off;
    
    %% 2D RGC RESPONSE
    figure(4)
    subplot(1,length(ratios),r); cla
    plot(filteredSubsampled(midpoint_rsp,:), 'k', 'LineWidth',4); colormap gray; 
    xlabel('x position (deg)')
    ylabel('modulation (a.u.)')
    title(labelsFilter)
    set(gca, 'XLim', [0 size(filteredSubsampled,1)], 'TickDir', 'out', 'FontSize', 15), ...
%         'XTick', [0,midpoint_rsp,size(filteredSubsampled,1)], ...
%         'XTickLabel', {'-1', '0', '1'});
    axis square; box off;
    
    %% 1D DOG FFT
    
    % Get 1D x axis
    ncones = length(G_1D);
    x      = (0:ncones-1)/ncones*rgcParams.fov; % degrees
    fs     = (0:ncones-1)/2;
    
    figure(5); clf; hold all;
    plot(fs-(max(fs)/2), G_1D, 'color', colors(r,:), 'LineWidth', 4);
    
    %% 2D DOG FFT
    figure(6)
    subplot(1,length(ratios),r); cla; hold on;
    imagesc(fftshift(G_2D)); colormap gray;
    
    % draw center cross
    plot([midpoint, midpoint], [1,size(G_2D,1)], ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);
    plot([1,size(G_2D,1)],[midpoint, midpoint],  ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);

    % draw stim sf cross
    plot(xvals, yvals, 'r:', 'lineWidth',3);
    colorbar; 
    set(gca, 'CLim', [0 0.7], 'XLim', [0 size(G_2D,1)])
    xlabel('cycles/deg');
    ylabel('cycles/deg');
    title(labelsFilter)
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
       'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'XTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)}, ...
        'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'YTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)});
    axis square; box off;
    
    %% 1D RGC RESPONSE
    figure(7); 
    
    ncones = length(filteredSubsampled_amps_1D);
    x      = (0:ncones-1)/ncones*rgcParams.fov; % degrees
    fs     = (0:ncones-1)/2;
    
    subplot(1,length(ratios),r); cla; hold on;
    plot(fs, fftshift(filteredSubsampled_amps_1D), 'k', 'LineWidth', 2); 
%     xlim([0 max(fs)/2])
    plot(rgcParams.stimSF * [1 1], [0 1], 'r:', 'LineWidth', 4);
    xlabel('sf (cycles/deg)')
    ylabel('modulation (a.u.)')
    title(labelsFilter)
%     set(gca, 'YLim', [0 2], ...
%         'TickDir', 'out', 'FontSize', 15, 'XTick', [0, midpoint], ...
%         'XTickLabel', {'0', num2str(midpoint/rgcParams.fov)});
%     axis square; box off;
%     
    %% 2D RGC RESPONSE FFT
    figure(8)
    subplot(1,length(ratios),r); cla
    filteredSubsampled_amps_2D_norm = fftshift(filteredSubsampled_amps_2D)./max(filteredSubsampled_amps_2D(:));
    imagesc(filteredSubsampled_amps_2D_norm); colormap gray;
    colorbar; 
    set(gca, ...
        'CLim', [-1 1])
    xlabel('sf (cycles/deg)')
    ylabel('sf (cycles/deg)')
    title(labelsFilter)
      
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [midpoint_rsp*0.5, midpoint_rsp, midpoint_rsp*1.5], ...
        'XTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)}, ...
        'YTick', [midpoint_rsp*0.5, midpoint_rsp, midpoint_rsp*1.5], ...
        'YTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)});
    axis square; box off;

end

figure(5);
 xlim([0 max(fs)/2])
hold on
plot(rgcParams.stimSF * [1 1], [0 .6], 'r:', 'LineWidth', 4);
xlabel('Spatial frequency (cycles/deg)'); 
ylabel('Normalized Modulation (a.u.)'); 
plot([0,sz(1)], [0 0], 'k'); hold on;
title(labelsFilter)
%     set(gca, 'YLim', [-0.05,1],'XLim', (midpoint + [-10, 10]), 'TickDir', 'out', 'FontSize', 15, 'XTick', [midpoint-10,midpoint,midpoint+10], ...
%         'XTickLabel', {'-0.25', '0', '0.25'});
axis square; box off;
        
figure(9);
plot([0,sz(1)], [0 0], 'k'); hold on;
xlabel('x position (deg)')
ylabel('modulation (a.u.)')
title('DoG filters')
set(gca,'YLim', [-0.05,1], 'XLim', [0, sz(1)], 'TickDir', 'out', 'FontSize', 15, 'XTick', [0, midpoint, sz(1)], ...
    'XTickLabel', {'-1', '0', '1'});
axis square; box off;
legend(labelsRatio, 'Location','Best'); legend boxoff;

hgexport(fH1, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs'))
hgexport(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'))
hgexport(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcResp1D'))
hgexport(fH4, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcResp2D'))
hgexport(fH5, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGFFT1D'))
hgexport(fH6, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGFFT2D'))
hgexport(fH7, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcRespFFT1D'))
hgexport(fH8, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcRespFFT2D'))
hgexport(fH9, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGAllFFT'))

print(fH1, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGs'), '-dpng')
print(fH2, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'Arrays'), '-dpng')
print(fH3, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcResp1D'), '-dpng')
print(fH4, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcResp2D'), '-dpng')
print(fH5, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGFFT1D'), '-dpng')
print(fH6, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGFFT2D'), '-dpng')
print(fH7, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcRespFFT1D'), '-dpng')
print(fH8, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'rgcRespFFT2D'), '-dpng')
print(fH9, fullfile(baseFolder, 'figures', 'RGC_layer', 'cartoons', 'DoGAllFFT'), '-dpng')



