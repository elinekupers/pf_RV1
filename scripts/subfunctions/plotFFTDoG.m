function [] = plotFFTDoG(DoGfilter,  rgcParams, expParams, saveFigs)



fftAmpsDoG = abs(fft2(DoGfilter));
fftAmpsDoG_normalized = fftAmpsDoG./max(fftAmpsDoG(:))';

midpoint = ceil(size(fftAmpsDoG,1)/2);
quarterpoint = midpoint/2;
sfOfStim = expParams.spatFreq*2;

% Visualize fft amps
fH = figure(1); clf;
imagesc(fftshift(fftAmpsDoG)); hold all;

% center cross
plot([midpoint, midpoint], [1,size(fftAmpsDoG,1)], 'w:');
plot([1,size(fftAmpsDoG,1)],[midpoint, midpoint],  'w:');

% stim sf cross
plot([midpoint-sfOfStim, midpoint-sfOfStim], [1,size(fftAmpsDoG,1)], 'r:', 'lineWidth',3);
plot([1,size(fftAmpsDoG,1)],[midpoint-sfOfStim, midpoint-sfOfStim],  'r:', 'lineWidth',3);

plot([midpoint+sfOfStim, midpoint+sfOfStim], [1,size(fftAmpsDoG,1)], 'r:', 'lineWidth',3);
plot([1,size(fftAmpsDoG,1)],[midpoint+sfOfStim, midpoint+sfOfStim],  'r:', 'lineWidth',3);

colormap gray; colorbar; axis square;

title(sprintf('Normalized FFT amps DoG filter %d',rgcParams.cone2RGCRatio), 'FontSize',17);
set(gca,'CLim', [0 max(fftAmpsDoG(:))], ...
    'FontSize', 17, 'TickDir', 'out', ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{num2str(-quarterpoint), '0', num2str(quarterpoint)}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{num2str(-quarterpoint), '0', num2str(quarterpoint)});
xlabel('cycles/deg'); box off;


if saveFigs
    figure(fH);
    hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('rgcDoGFilter_ratio%d', rgcParams.cone2RGCRatio)))
end

return