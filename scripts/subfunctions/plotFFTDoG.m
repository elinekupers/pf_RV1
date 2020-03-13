function [] = plotFFTDoG(DoGfilter,  rgcParams, expParams, saveFigs)

fftAmpsDoG = abs(fft2(DoGfilter));
fftAmpsDoG_normalized = fftAmpsDoG./max(fftAmpsDoG(:))';

midpoint = ceil(size(fftAmpsDoG,1)/2);
quarterpoint = midpoint/2;

% get stim sf and get coords to draw circle
sfOfStim = expParams.spatFreq * rgcParams.fov;

th = 0:pi/50:2*pi;
xvals = sfOfStim * cos(th) + midpoint;
yvals = sfOfStim * sin(th) + midpoint;


% Visualize fft amps
fH = figure(1); clf;
imagesc(fftshift(fftAmpsDoG_normalized)); hold all;

% draw center cross
plot([midpoint, midpoint], [1,size(fftAmpsDoG,1)], 'w:');
plot([1,size(fftAmpsDoG,1)],[midpoint, midpoint],  'w:');

% draw stim sf cross
plot(xvals, yvals, 'r:', 'lineWidth',3);

% Make figure pretty
colormap gray; colorbar; axis square;
title(sprintf('Normalized FFT amps DoG filter %d',rgcParams.cone2RGCRatio), 'FontSize',17);
set(gca,'CLim', [0 max(fftAmpsDoG_normalized(:))], ...
    'FontSize', 17, 'TickDir', 'out', ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{num2str(-quarterpoint/rgcParams.fov), '0', num2str(quarterpoint/rgcParams.fov)});
xlabel('cycles/deg'); box off;


if saveFigs
    figure(fH);
    hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('rgcDoGFilter_ratio%d', rgcParams.cone2RGCRatio)))
end

return