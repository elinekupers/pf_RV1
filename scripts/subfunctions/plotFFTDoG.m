function [] = plotFFTDoG(DoGfilter,  rgcParams)


% take fft, get amplitudes
fftAmpsDoG = abs(fft2(DoGfilter./sum(abs(DoGfilter(:)))));

midpoint = ceil(size(fftAmpsDoG,1)/2);
quarterpoint = midpoint/2;

% get stim sf and get coords to draw circle
sfOfStim = rgcParams.stimSF * rgcParams.fov;

th = 0:pi/50:2*pi;
xvals = sfOfStim * cos(th) + midpoint;
yvals = sfOfStim * sin(th) + midpoint;


% Visualize fft amps
fH = figure(1); clf; set(gcf, 'Position', [397, 748, 1163, 590], 'Color', 'w');

subplot(121)
imagesc(fftshift(fftAmpsDoG)); hold all;

% draw center cross
plot([midpoint, midpoint], [1,size(fftAmpsDoG,1)], ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);
plot([1,size(fftAmpsDoG,1)],[midpoint, midpoint],  ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth',4);

% draw stim sf cross
plot(xvals, yvals, 'r:', 'lineWidth',3);

% Make figure pretty
colormap gray; colorbar; axis square;
title(sprintf('FFT amps DoG filter %d',rgcParams.cone2RGCRatio), 'FontSize',17);
set(gca,'CLim', [0 max(fftAmpsDoG(:))], ...
    'FontSize', 17, 'TickDir', 'out', ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)});
    xlabel('cycles/deg'); box off;

    
% Get 1D representation
ncones = size(fftAmpsDoG,1);
x      = (0:ncones-1)/ncones*rgcParams.fov; % degrees
fs     = (0:ncones-1)/2;
G      = fftshift(fftAmpsDoG);

subplot(122); cla;
plot(fs-(max(fs)/2), G(:,midpoint)', 'k', 'LineWidth', 4); xlim([0 max(fs)/2])
hold on
plot(rgcParams.stimSF * [1 1], [0 1], 'r:', 'LineWidth', 4);
xlabel('Spatial frequency (cycles/deg)'); 
ylabel('Normalized Modulation (a.u.)'); 
set(gca, 'TickDir','out', 'FontSize', 17); box off;

if rgcParams.saveFigs
    figure(fH);
    hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('rgcDoGFilter_ratio%d', rgcParams.cone2RGCRatio)))
end

return