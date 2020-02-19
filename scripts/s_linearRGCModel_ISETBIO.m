%% s_linearRGCModel_ISETBIO
%
% First attempt script to build an RGC layer on top of the current
% computational observer model. The RGC layer is based on the retina-V1
% model published by Bradley, Abrams and Geisler (2014) in JoV.
% 
% Our RGC layer of the model has two stages, after getting the cone
% currents, first one is linear and the second is non-linear:
%   0.  Obtain cone current, either by running observer model, or loading
%   in saved cone currents. 
%   1.  Cone current -> filtered by Difference of Gaussians (DoG) 
%   2.  Filtered response --> sub sampling
%
%
% To recompute cone currents without eyemovements or Gabor stim phase
% shifts, use the following command:
% runComputationalObserverModel('defaultnophaseshift', 'saveFolder', ...
%                             'conecurrentRV1','seed',1,'currentFlag',true)

%% 0. Load cone current

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/'; % can also be ogRootPath if rerunning model

% Get experimental params
subFolder = 'conecurrentRV1';
expName   = 'defaultnophaseshift';
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)

fName     = 'current_OGconeOutputs_contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat';
load(fullfile(baseFolder, 'data', expName, subFolder, fName));

% Take one trial, one orientation, one phase, and mean across time points
img      = squeeze(mean(current(1,:,:,:,1),4));

% Plot cone current retinal image
figure(1); clf; imagesc(img); 
axis square; colormap gray; set(gca, 'CLim', [min(img(:)),0], 'TickDir', 'out')
title('Cone current, 1 trial, average across time points'); colorbar;
xlabel('X cones'); ylabel('Y cones'); 


%% OBSOLETE: Get cone density and spacing of photoreceptors

 % Convertion deg to m
deg2m    = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter

% Compute x,y position in m of center of retinal patch from ecc and angle
whichEye = 'left'; 
[x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
x = x * deg2m;  y = y * deg2m;
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

% Set the field of view (degrees)
cMosaic.setSizeToFOV(cparams.cmFOV);

% Get cone density from eccentricity stored in cMosaic object.
coneDensityCountsPerDeg2 = eccen2density(cMosaic, 'deg');

% Convert cone density (counts/deg2) to cone spacing (count/deg
spacingCountsPerDeg = sqrt(coneDensityCountsPerDeg2)/coneDensityCountsPerDeg2;

%% Get cone to RGC ratio

% OBSOLETE: RGC ratio from Watson (2014) model
% WatsonRGCCalc       = WatsonRGCModel('generateAllFigures', false);
% [xPos, yPos]        = pol2cart(cparams.polarAngle, cparams.eccentricity);
% eccXYDegs           = [xPos(:) yPos(:); 0 0];
% mRGCRFtoConesRatios = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYDegs, 'left');
% mRGCspacing         = (1/mRGCRFtoConesRatios(1)) * spacingCountsPerDeg;

% Take several cone 2 RGC density ratio's

cone2RGCRatios = [1,2,4,8,9]; %[1,2,4,8,9,16,18];

% Nr of cones in mosaic
sz    = size(img); 
cx    = sz(1);      % #cones on x-axis
cy    = sz(2);      % #cones on y-axis

[X,Y] = meshgrid(1:cx,1:cy);

% DoG Params 
kc = 1;    % Gauss center sigma. (as in Bradley et al. 2014 paper)
ks = 2;    % Gauss surround sigma. Range: ks > kc. Note: Bradley paper uses 10.1, seems very large to us 

wc = 0.53; % DoG center Gauss weight. Range: [0,1]. Currently not used
ws = 1-wc; % DoG surround Gauss weight.  Currently not used

% Preallocate space
filteredConeCurrent = NaN(length(cone2RGCRatios), cx, cy);
rgcResampled        = cell(size(cone2RGCRatios));

% Set subplots
figure(1); clf; set(gcf, 'Position', [244,46,2316,1299], 'Color', 'w');
cols    = length(cone2RGCRatios);
rows    = 3;


for ii = 1:length(cone2RGCRatios)
    
    % Center Gauss RGC
    sigma.center = cone2RGCRatios(ii);
    
    % Surround Gauss RGC
    sigma.surround = sigma.center*ks;
    
    % ratio center surround
    sigma.ratio = sigma.surround/sigma.center;
    
    % ratio of the surround volume to the center volume
    vol.ratio = 1;
    
    % Create RGC grid by resampling with cone2RGC ratio.
    rgc = getRGCMosaic(cx, cy, cone2RGCRatios(ii));
    rgcMask = logical(rgc);
    
    [f,xx,yy] = makedog2d(cx,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);
    
    % Convolve image with DoG filter
    filteredConeCurrent(ii,:,:) = conv2(img, f, 'same');
    
    % Resample RGC image
    tmpIm = squeeze(filteredConeCurrent(ii, :, :));
    rgcIm = rgc;
    rgcIm(rgcMask) = tmpIm(rgcMask);
    rgcResampled{ii} = reshape(rgcIm, cx, cy);
    
    % Plot mosaics
    figure(1); hold all;
    subplot(rows, cols, ii); 
    plot(X, Y, 'ko'); hold on;
    plot(X(logical(rgc)), Y(logical(rgc)), 'r.'); 
    xlim([20 60]); ylim([20 60]); axis square; box off;
    title(sprintf('Sub sampling, Cones:RGCs=1:%d', cone2RGCRatios(ii)));
    xlabel('# cones');
    ylabel('# cones');
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(rows, cols, (ii+cols));
    %imagesc(f, [-1 1]);  colorbar;
    mesh(xx,yy,f); axis square; box off; view([0 0]);
    title(sprintf('DoG filter 1:%d', cone2RGCRatios(ii)));
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('normalized response')
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    subplot(rows, cols, (ii+(cols*2)));
    imagesc(rgcResampled{ii}); colormap gray; colorbar; axis square; box off;
    title(sprintf('Filtered cone current (pA) 1:%d', cone2RGCRatios(ii)));
    xlabel('# cones (x)');
    ylabel('# cones (y)');
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
end

%% Get cone to retinal ganglion cell ratio to do subsampling



