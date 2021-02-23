load('/Volumes/server/Projects/PerformanceFields_RetinaV1Model/data/default/conecurrentRV1/OGconeOutputs_contrast1.0000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat')
coneresp = squeeze(absorptions(1,:,:,1,1));

rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;               % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight.
rgcParams.cRows      = size(coneresp,1);
rgcParams.cCols      = size(coneresp,2);
rgcParams.seed       = 1;

figure;
for ii = 1:5;
    
rgcParams.cone2RGCRatio = ii;
% Center Gauss RGC
sigma.center = rgcParams.cone2RGCRatio*rgcParams.DoG.kc;

% Surround Gauss RGC
sigma.surround = rgcParams.cone2RGCRatio*rgcParams.DoG.ks;

% ratio center surround
sigma.ratio = sigma.surround/sigma.center;

% ratio of the surround volume to the center volume
vol.ratio = rgcParams.DoG.ws/rgcParams.DoG.wc;

% Create RGC grid by resampling with cone2RGC ratio.
rowIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cRows;
colIndices = 1:rgcParams.cone2RGCRatio:rgcParams.cCols;
conearray = zeros(rgcParams.cCols, rgcParams.cRows);

rgcarray = conearray;
rgcarray(rowIndices, colIndices) = 1;

% Create DoG filter
[DoGfilter,xx,yy] = makedog2d(rgcParams.cRows,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);

% mean cone absorption--> get poission noise to pad surround
mnPad = ones(size(coneresp,1)*2,size(coneresp,2)*2).*mean(coneresp(:));
mnPadPoisson = iePoisson(mnPad, 'noiseFlag', 'random', 'seed', rgcParams.seed);
rowStart = floor(size(coneresp,1)/2)+1;
colStart = floor(size(coneresp,2)/2)+1;
imgMnPadded = mnPadPoisson;
imgMnPadded([rowStart:(rowStart+size(coneresp,1)-1)],[colStart:(colStart+size(coneresp,2)-1)]) = coneresp;

% Convolve padded image with DoG filter
filteredConeRespFull = conv2(imgMnPadded, DoGfilter, 'same');
filteredConeRespPadded = filteredConeRespFull([rowStart:(rowStart+size(coneresp,1)-1)],[colStart:(colStart+size(coneresp,2)-1)]);

% Convolve default image with DoG filter
filteredConeResp = conv2(coneresp, DoGfilter, 'same');

% Resample RGC images
rgcResponsePadded = squeeze(filteredConeRespPadded(rowIndices, colIndices));
rgcResponse = squeeze(filteredConeResp(rowIndices, colIndices));

 % we divide by 2, because the time sampling is at 2 ms.
subplot(2,5,ii);
imagesc(rgcResponsePadded./2); colormap gray; axis square; box off;
title(sprintf('Mean RGC response using padded image, ratio %d', ii));
xlabel('# cones (rows)');
ylabel('# cones (cols)');

subplot(2,5,ii+5);
imagesc(rgcResponse./2); colormap gray; axis square; box off;
title(sprintf('Mean RGC response no padding, ratio %d', ii));
xlabel('# cones (rows)');
ylabel('# cones (cols)');
end