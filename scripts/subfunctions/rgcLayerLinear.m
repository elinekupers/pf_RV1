function [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(coneData, rgcParams)
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
%   2.  Filtered response --> sub sampling (Non-linear)
%
% INPUTS:
%   coneData            : 5 dimensional array of cone currents 
%                           (trials x rows x cols x time x stim phase)
%   rgcParams           : params that define the DoG filter and the
%                           cone:RGC ratio
% 
% OUTPUTS:
%   rgcParams           : RGC
%   rgcarray            :
%   DoGfilter           :
%   filteredConeCurrent :


% reshape cone data 
% original: trials x rows x cols x time x stim phase
% reshaped: rows x cols x time x all trials    
permutedConeData = permute(coneData, [2, 3, 4, 1, 5]);
reshapedConeData = reshape(permutedConeData, rgcParams.cRows, rgcParams.cCols, rgcParams.timePoints, []);

% Create cone array
conearray = zeros(rgcParams.cCols, rgcParams.cRows);

% Center Gauss RGC
sigma.center = rgcParams.cone2RGCRatio;

% Surround Gauss RGC
sigma.surround = sigma.center*rgcParams.DoG.ks;

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
        
        img =  reshapedConeData(:,:,t, ii);

        % Convolve image with DoG filter
        filteredConeCurrent = conv2(img, DoGfilter, 'same');

        % Resample RGC image
        rgcResponse(:,:,t,ii) = squeeze(filteredConeCurrent(rowIndices, colIndices));
    end
end

% reshape rgc response back to original dimensions
[numRGCRows, numRGCCols, numTimePoints, ~]  = size(rgcResponse);
rgcResponse = reshape(rgcResponse, numRGCRows, numRGCCols, numTimePoints, rgcParams.nTrials, []);
rgcResponse = permute(rgcResponse, [4, 1, 2, 3, 5]);


if rgcParams.verbose
    [X,Y] = meshgrid(1:rgcParams.cRows,1:rgcParams.cCols);
    rgcResponse_mn = squeeze(mean(mean(rgcResponse(:,:,:,:,1),1),4));
    fH = figure(99); clf; set(gcf, 'Position', [244,680,2316,665], 'Color', 'w');
    
    % Plot mosaics
    subplot(1, 3, 1); hold all;
    plot(X,Y, 'k.', 'MarkerSize',9); hold on;
    plot(X(logical(rgcarray)),Y(logical(rgcarray)), 'r.');
    axis square
    title(sprintf('Sub sampling, cones:RGCs=1:%d',rgcParams.cone2RGCRatio));
    xlabel('# cones (rows)');
    ylabel('# cones (cols)');
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(1, 3, 2);
    imagesc(DoGfilter, [-1 1]);  colorbar;
    surf(xx,yy,DoGfilter); axis square; box off; view([0 0]);
    title(sprintf('DoG filter, center:surround=%d:%2.2f',sigma.center,sigma.surround));
    xlabel('x-axis (# cones)');
    ylabel('y-axis (# cones)');
    zlabel('normalized response')
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(1, 3, 3);
    imagesc(rgcResponse_mn); colormap gray; axis square; box off;
    title('Mean RGC response across trials and time');
    xlabel('# cones (rows)');
    ylabel('# cones (cols)');
    h = colorbar; ylabel(h, 'current (pA)')
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    if rgcParams.saveFigs
        hgexport(fH, fullfile(pfRV1rootPath, 'figures', sprintf('rgclayerresponse_%d',rgcParams.cone2RGCRatio)));
    end
end

return