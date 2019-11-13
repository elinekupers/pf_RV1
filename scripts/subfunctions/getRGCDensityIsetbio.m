function plotRGCDensityIsetbio(fov)
%% NOT YET IMPLEMENTED: plotConeDensityIsetbio(angDeg, eccMM, sourceDataset)
% 
% function to plot mRGC density using ISETBIO toolbox
% 
% INPUTS:
%   fov           :   integer to define field of view of RGC mosaic (in degrees)


% Check inputs
if ~isempty(fov, 'var') || ~exist('fov','var')
    fov = 2;
end

% We first need to create a cone mosaic object
cMosaic = coneMosaic('center', [0 0]);
cMosaic.setSizeToFOV(fov);

% Create bipolar layer
bpL = bipolarLayer(cMosaic);

bpMosaicParams.spread  = 2;  % ???RF diameter w.r.t. input samples 
bpMosaicParams.stride  = 2;  % ???RF diameter w.r.t. input samples
bpL.mosaic{1} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);

% Create a rgc layer
rgcL = rgcLayer(bpL);

% Define rfDiameter? (Spread and stride are not working?)
rgcParams.rfDiameter = 2;

% Get mosaic
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget',rgcParams);

% Get nr of cells
[x_cells, y_cells, ~] = size(rgcL.mosaic{1}.cellLocation);
nrOfCells = prod([x_cells,y_cells]);

% What's the size of the rgc mosaic patch?
patchSize  = prod(rgcL.size);

% Density within fov
rgcDensity = nrOfCells * (1/patchSize);


