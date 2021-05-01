function [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(coneData, rgcParams, expParams)
% Function that represents a simple RGC layer with linear-nonlinear cascade:
% RGC layer takes as input the cone current from our previous ISETBIO
% computational observer model.
%
% [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = ...
%                               rgcLayerLinear(coneData, rgcParams)
%
% The RGC layer is based on the retina-V1 model published by Bradley,
% Abrams and Geisler (2014) in JoV, and has two stages:
%
%   1.  Cone current -> filtered by Difference of Gaussians (DoG) (Linear)
%   2.  Filtered response --> sub sampling (Linear)
%
% INPUTS:
%   coneData            : 5 dimensional array of cone currents
%                           (trials x rows x cols x time x stim phase)
%   rgcParams           : params that define the DoG filter and the
%                           cone:RGC ratio
%   expParams           : structure created by loadExpParams to describe
%                           experiment
%
% OUTPUTS:
%   rgcResponse         : mRGC responses (i.e. filtered and resampled cone
%                           data), 5-D array
%   rgcarray            : Resampling matrix of mRGC layer (same size as cone
%                           matrix)
%   DoGfilter           : Difference of Gaussian filter used to filter cone
%                           data
%   filteredConeCurrent : Filtered, but not resampled, mRGC data.
%
%
% Written by EK @ NYU, 2020

%% reshape cone data
% original: trials x rows x cols x time x stim phase
% reshaped: rows x cols x time x all trials
permutedConeData = permute(coneData, [2, 3, 4, 1, 5]);
reshapedConeData = reshape(permutedConeData, rgcParams.cRows, rgcParams.cCols, rgcParams.timePoints, []);

% Create cone array
conearray = zeros(rgcParams.cCols, rgcParams.cRows);

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

rgcarray = conearray;
rgcarray(rowIndices, colIndices) = 1;

% Create DoG filter
[DoGfilter,xx,yy] = makedog2d(rgcParams.cRows,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);

for ii = 1:size(reshapedConeData,4)
    
    for t = 1:size(reshapedConeData,3)
        
        % get single time frame image
        img =  reshapedConeData(:,:,t, ii);
        
        % mean cone absorption--> get poission noise to pad surround
        mnPad = ones(size(img,1)*2,size(img,2)*2).*mean(img(:));
        if strcmpi(expParams.cparams.noise, 'none')
            mnPadPoisson = mnPad;
        else
            mnPadPoisson = iePoisson(mnPad, 'noiseFlag', 'random', 'seed', rgcParams.seed);
        end
        rowStart = floor(size(img,1)/2)+1;
        colStart = floor(size(img,2)/2)+1;
        imgMnPadded = mnPadPoisson;
        imgMnPadded([rowStart:(rowStart+size(img,1)-1)],[colStart:(colStart+size(img,2)-1)]) = img;
        
        % Convolve image with DoG filter
        filteredConeCurrentFull = conv2(imgMnPadded, DoGfilter, 'same');
        filteredConeCurrent = filteredConeCurrentFull([rowStart:(rowStart+size(img,1)-1)],[colStart:(colStart+size(img,2)-1)]);
        
        % Resample RGC image
        rgcResponse(:,:,t,ii) = squeeze(filteredConeCurrent(rowIndices, colIndices));
    end
end

% reshape rgc response back to original dimensions
[numRGCRows, numRGCCols, numTimePoints, ~]  = size(rgcResponse);
rgcResponse = reshape(rgcResponse, numRGCRows, numRGCCols, numTimePoints, rgcParams.nTrials, []);
rgcResponse = permute(rgcResponse, [4, 1, 2, 3, 5]);


%% Some visualization, if verbose = true
if rgcParams.verbose
    
    % Plot fft amps of filter
    plotFFTDoG(DoGfilter,  rgcParams)
    
    % Plot array, filter and response
    [X,Y] = meshgrid(1:rgcParams.cRows,1:rgcParams.cCols);
    rgcResponse_mn = squeeze(mean(mean(rgcResponse(:,:,:,:,1),1),4));
    
    fH = figure(99); clf; set(gcf, 'Position', [244,680,2316,665], 'Color', 'w');
    
    % Plot mosaics
    subplot(1, 4, 1); hold all;
    plot(X,Y, 'k.', 'MarkerSize',9); hold on;
    plot(X(logical(rgcarray)),Y(logical(rgcarray)), 'r.');
    axis square
    title(sprintf('Sub sampling, cones:RGCs=1:%1.1f',rgcParams.cone2RGCRatio));
    xlabel('# cones (rows)');
    ylabel('# cones (cols)');
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(1, 4, 2);
    imagesc(DoGfilter, [-0.2, 1]);  colorbar; colormap gray;
    title(sprintf('Top view DoG, size center:surround=%2.2f:%2.2f',sigma.center,sigma.surround));
    xlabel('rows (# cones)'); ylabel('cols (# cones)');
    set(gca, 'TickDir', 'out', 'FontSize', 12); axis square; box off;

    subplot(1, 4, 3);
    surf(xx,yy,DoGfilter); axis square; box off; view([0 0]);
    
    title(sprintf('Side view DoG, size center:surround=%2.2f:%2.2f',sigma.center,sigma.surround));
    xlabel('x-axis (# cones)');
    ylabel('y-axis (# cones)');
    zlabel('normalized response')
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(1, 4, 4); % we divide by 2, because the time sampling is at 2 ms.
    imagesc(rgcResponse_mn./2); colormap gray; axis square; box off;
    title('Mean RGC response across trials and time');
    xlabel('# cones (rows)');
    ylabel('# cones (cols)');
    if strcmp(rgcParams.inputType, 'current')
        yl = 'current (pA)';
    else
        yl = 'absorption rate (photons/ms)';
    end
    h = colorbar; ylabel(h, yl);
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    if rgcParams.saveFigs
        if ~exist(fullfile(pfRV1rootPath, 'figures', sprintf('RGC_%s', rgcParams.subFolder))),
            mkdir(fullfile(pfRV1rootPath, 'figures', sprintf('RGC_%s', rgcParams.subFolder)));
        end
        hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('RGC_%s', rgcParams.subFolder), sprintf('rgclayerresponse_%d_contrast%1.4f_%s.eps',rgcParams.cone2RGCRatio, rgcParams.c, rgcParams.inputType)));
    end
    
end

return