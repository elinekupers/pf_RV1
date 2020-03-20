function [] = plot1DFFTDoG()

% Get RGC layer params
% Take several cone:RGC density ratio's (for subsampling)
cone2RGCRatios = 1:1:5; % linear ratio

doSubsampling  = true;

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = true; % print figures or not
rgcParams.saveFigs   = true;

% Define DoG Params
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = 10.1/3;             % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10)
rgcParams.DoG.wc     = 0.53;               % DoG center Gauss weight. Range: [0,1].
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.cRows      = 79;  % #cones on y-axis
rgcParams.cCols      = 79;  % #cones on x-axis
rgcParams.fov        = 2;   % deg
rgcParams.stimSF     = 4;   % cpd

% get stim sf and get coords to draw circle
sfOfStim = rgcParams.stimSF * rgcParams.fov;
th = 0:pi/50:2*pi;

cmap = parula(length(cone2RGCRatios)+1);

% Get figure
fH = figure(1); clf; set(gcf, 'Position', [228,218,721,597], 'Color', 'w'); hold all;

for ii = 1:length(cone2RGCRatios)
    
    % set cone:RGC ratio
    rgcParams.cone2RGCRatio = cone2RGCRatios(ii);
    
    % Center Gauss RGC
    sigma.center = rgcParams.cone2RGCRatio*rgcParams.DoG.kc;
    
    % Surround Gauss RGC
    sigma.surround = rgcParams.cone2RGCRatio*rgcParams.DoG.ks;
    
    % ratio center surround
    sigma.ratio = sigma.surround/sigma.center;
    
    % ratio of the surround volume to the center volume
    vol.ratio = rgcParams.DoG.ws/rgcParams.DoG.wc;
    
    % Create cone array
    conearray = zeros(rgcParams.cCols, rgcParams.cRows);
    
    % Create RGC grid by resampling with cone2RGC ratio.
    rowIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cRows;
    colIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cCols;
    
    % Create RGC array
    rgcarray = conearray;
    
    if doSubsampling
        rgcarray(rowIndices, colIndices) = 1;
    end
    
    % Create DoG filter
    [DoGfilter,xx,yy] = makedog2d(rgcParams.cRows,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);
    
    % normalize DoG, take fft, get amplitudes
    fftAmpsDoG = abs(fft2(DoGfilter./sum(abs(DoGfilter(:)))));
    
    midpoint = ceil(size(fftAmpsDoG,1)/2);
    
    % Visualize  1D representation
    ncones = size(fftAmpsDoG,1);
    fs     = (0:ncones-1)/2;
    G      = fftshift(fftAmpsDoG);
    
    plot(fs-(max(fs)/2), G(:,midpoint)', 'color', cmap(ii,:,:), 'LineWidth', 4);
    labels{ii} = sprintf('Cone:RGC = 1:%1.1f', rgcParams.cone2RGCRatio);
end

plot(rgcParams.stimSF * [1 1], [0 1], 'k:', 'LineWidth', 4);
xlabel('Spatial frequency (cycles/deg)');
ylabel('Normalized Modulation (a.u.)');
title('1D representation of DoG RGC filter')
set(gca, 'TickDir','out', 'FontSize', 17); box off;
xlim([0 max(fs)/2]); axis square
legend(labels, 'Location', 'EastOutside'); legend boxoff;

if rgcParams.saveFigs
    figure(fH);
    print(fH,fullfile(pfRV1rootPath, 'figures', sprintf('rgcDoGFilter_1D_ratio1-5_subsampling%d',doSubsampling)), '-dpng')
    hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('rgcDoGFilter_1D_ratio1-5_subsampling%d', doSubsampling)))
end
