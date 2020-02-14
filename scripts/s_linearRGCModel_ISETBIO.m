%% s_linearRGCModel_ISETBIO


% load cone current params
subFolder = 'conecurrentRV1';
expName   = 'defaultnophaseshift';
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)

fName     = 'current_OGconeOutputs_contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat';
load(fullfile(ogRootPath, 'data', expName, subFolder, fName));

% take one trial, one orientation, one phase, and mean across time points
img      = squeeze(mean(current(1,:,:,:,1),4));

figure; imagesc(img); 
axis square; colormap gray; set(gca, 'CLim', [min(img(:)),0], 'TickDir', 'out')
title('Cone current, 1 trial, average across time points'); colorbar;
xlabel('X cones'); ylabel('Y cones'); 


%% Get cone density and spacing of photoreceptors

 % Convertion deg to m
deg2m    = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
whichEye = 'left'; 
% Compute x,y position in m of center of retinal patch from ecc and angle
[x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
x = x * deg2m;  y = y * deg2m;
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

% Set the field of view (degrees)
cMosaic.setSizeToFOV(cparams.cmFOV);
% propCovered = getBanks1991ConeCoverage(cparams.eccentricity); % proportion
% cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered; % meters
% cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered; % meters

coneDensityCountsPerDeg2 = eccen2density(cMosaic, 'deg');

spacingCountsPerDeg = sqrt(coneDensityCountsPerDeg2)/coneDensityCountsPerDeg2;

%% Get cone to RGC ratio

WatsonRGCCalc       = WatsonRGCModel('generateAllFigures', false);
[xPos, yPos]        = pol2cart(cparams.polarAngle, cparams.eccentricity);
eccXYDegs           = [xPos(:) yPos(:); 0 0];
mRGCRFtoConesRatios = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYDegs, 'left');

mRGCspacing = (1/mRGCRFtoConesRatios(1)) * spacingCountsPerDeg;

%% Center Gauss RGC

sigma.center = mRGCspacing;

%Optics of the human eye
sz    = size(img); 
degX  = sz(1)/ppd;
degY  = sz(2)/ppd;

[Y,X] = meshgrid(1:sz(2),1:sz(1));

% Get center Gaussian of specified standard deviation
GausCenter = fft_Gaus(g_stds,sqrt(((X-center(sz(1))).^2)./(degX.^2)+((Y-center(sz(2))).^2)./(degY.^2)));

% Convolve with image in fourier space, then transform back
centerGausImg   = fftshift(ifft2(fft_img.*ifftshift(GausCenter)));

%% Surround Gauss RGC 
% computation from bradley et al. 2014 paper

kc = 1;    % Gauss center sigma.
ks = 10.1; % Gauss surround sigma. Range: ks > kc
wc = 0.53; % DoG center Gauss weight. Range: [0,1]
ws = 1-wc; % DoG surround Gauss weight

% The kth level of the output stack will estimate the result of convolving
% I with a Gaussian whose standard deviation is
%(log2(interp_scalar*2^k))/120 degrees.
k = 1; % ??? stack level?

% Surround SD
sdSurround = (log2(ks*2^k))/ppd;

% Gauss Surround
GausSurround = fft_Gaus(sdSurround,sqrt(((X-center(sz(1))).^2)./(degX.^2)+((Y-center(sz(2))).^2)./(degY.^2)));

% Again convolve with image in fourier space, then transform back
surroundGausImg   = fftshift(ifft2(fft_img.*ifftshift(Surround)));

% Take diff of Gaussians
imgDoG = P.wc.*T - (1-P.wc).*T_sur;


%% Get cone to retinal ganglion cell ratio to do subsampling



