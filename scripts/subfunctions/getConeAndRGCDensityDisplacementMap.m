function [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
            getConeAndRGCDensityDisplacementMap(sampleResEccen, sampleResPolAng, maxEccen, ...
            pixelsPerDegVisual, saveFigures, saveData)
%% getConeAndRGCDensityDisplacementMap(angDeg, eccMM, sourceDataset)
% 
% function to plot cone density using ISETBIO toolbox
% 
% INPUTS:
%   sampleResEccen      :   sampling resolution of eccentricities (in deg)
%   sampleResPolAng     :   sampling resolution of polar angles (in deg)
%   maxEccen            :   max eccentricity to plot (in deg)
%   pixelsPerDegVisual  :   number of pixels per deg visual angle, for
%                               plotting
%   saveFigures         :   save figures or not (bool), default: false
%   saveData            :   save data or not (bool), default: false
%
% OUTPUTS:
%   allMaps             :   cell with cone density map {1},
%                                   rgc density map {2},
%                                   midget RGC RF {3}, 
%                                   ratio midget RGC RF to cone density {4}



% Check inputs
if isempty(sampleResEccen) || ~exist('sampleResEccen','var')
    sampleResEccen = 0.01; % Needs to be 0.01, otherwise there is a chance of introducing an 'out of bound optimization error'
end

if isempty(sampleResPolAng) || ~exist('sampleResPolAng','var')
    sampleResPolAng = 5;
end

if isempty(maxEccen) || ~exist('maxEccen','var')
    maxEccen = 40; % Needs to be at least 40, otherwise there is a chance of introducing an 'out of bound optimization error'
end

if isempty(pixelsPerDegVisual) || ~exist('pixelsPerDegVisual','var')
    pixelsPerDegVisual = 10;
end

if isempty(saveFigures) || ~exist('saveFigures','var')
    saveFigures = false;
end

if isempty(saveData) || ~exist('saveData','var')
    saveData = false;
end

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end

dataDir = fullfile(pfRV1rootPath, 'external', 'data');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end



%% LOAD CONE AND RGC DATA USING RGC DISPLACEMENT MAP

% Define the cone and rgc density data to operate upon
subjectName = 'computedAverage';
coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',subjectName,'.mat']);
rgcDensityDataFileName  = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',subjectName,'.mat']);

% Create the displacement model
[ rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegVisual, opticDiscLocationByMeridian, mRF_RingCumulativeByMeridian, mRGC_RingCumulativeByMeridian, fitParamsByMeridian, ~, ~ ] = ...
    createDisplacementModel('verbose', true, ...
    'sampleResolutionDegVisual', sampleResEccen, ...
    'maxModeledEccentricityDegVisual', maxEccen, ...
    'meridianAngleResolutionDeg', sampleResPolAng, ...
    'displacementMapPixelsPerDegVisual', pixelsPerDegVisual);


% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngleSupport)
    fprintf('.');
    % obtain cone density
    fitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(meridianAngleSupport(mm));
    coneDensitySqDeg = zeroOpticDiscPoints(fitConeDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    coneDensityByMeridian(mm,:) = coneDensitySqDeg;
        
    % obtain the RGC density
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(meridianAngleSupport(mm));
    rgcDensitySqDeg = zeroOpticDiscPoints(fitRGCDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    rgcDensityByMeridian(mm,:) = rgcDensitySqDeg;
    
    % obtain the mRF density
    [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
    mRFDensityByMeridian(mm,:) = mRFDensitySqDeg;
    mRFtoConeDensityByMeridian(mm,:) = mRFtoConeDensityRatio;
    
end
fprintf('\n');

%% DEFINE OPTIC DISK AND SAMPLE BASE 

% Define the image sample base for transform from polar coords
imRdim = (max(regularSupportPosDegVisual) * pixelsPerDegVisual * 2) -1;

% Show and save the maps
polarMapNameList = {...
    'coneDensityByMeridian',...
    'rgcDensityByMeridian',...
    'mRFDensityByMeridian',...
    'mRFtoConeDensityByMeridian'};

% Obtain the boundary of the optic disc in image space. We have to reverse
% the boundary array for the y-axis to match our image convention
opticDiscLocationsImage     = convertPolarMapToImageMap(opticDiscLocationByMeridian, 'imRdim', imRdim);
opticDiscBoundary           = bwboundaries(imbinarize(opticDiscLocationsImage),'noholes');
opticDiscBoundaryArray      = opticDiscBoundary{1};
opticDiscBoundaryArray(:,1) = imRdim(1)-opticDiscBoundaryArray(:,1);
dashIndices                 = 1:10:size(opticDiscBoundaryArray,1);

% preallocate space
allMaps = cell(1,length(polarMapNameList));

%% VISUALIZE DATA

% loop over the maps
for vv = 1:length(polarMapNameList)
    
    % convert polar map to image
    mapImage = feval('convertPolarMapToImageMap', eval(polarMapNameList{vv}), 'imRdim', imRdim);
    
    % store data and sample resolutions
    allMaps{vv}.data = mapImage;
    allMaps{vv}.name = polarMapNameList{vv};
    allMaps{vv}.regularSupportPosDegVisual = regularSupportPosDegVisual;
    allMaps{vv}.sampleResPolAng = 0:sampleResPolAng:(360-1);
    allMaps{vv}.pixelsPerDegVisual = pixelsPerDegVisual;
    allMaps{vv}.maxEccen = maxEccen;
    
    % plot it
    fH = figure();
    fH.Renderer='Painters';
    climVals = [0,ceil(max(max(mapImage)))];
    tmp = strsplit(polarMapNameList{vv},'ByMeridian');
    titleString=tmp{1};
    
    % Handle the special case of the clim for mRFtoConeDensityByMeridian
    if strcmp(titleString,'mRFtoConeDensity')
        climVals(2)=2;
    end
    
    displayRetinalImage(mapImage, climVals, pixelsPerDegVisual, maxEccen, titleString);
    title(titleString,'FontSize',20);
    hold on
    
    % Plot optic disc boundaries with no data
    if median(arrayfun(@(x) mapImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
        opticDiscBorderColor = [0.25 0.25 0.25];
    else
        opticDiscBorderColor = [0.75 0.75 0.75];
        
    end
    plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
    
    if saveFigures
        % Save matlab fig and pdf
        savefig(fH, fullfile(figureDir, polarMapNameList{vv}))
        print(fullfile(figureDir, polarMapNameList{vv}), '-dpdf', '-fillpage')
    end
end

% Save data
if saveData
    save(fullfile(dataDir', 'coneDensityByMeridian.mat'), 'coneDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng');
    save(fullfile(dataDir, 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng');
    save(fullfile(dataDir', 'mRFtoConeDensityByMeridian.mat'),'mRFtoConeDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng');
    save(fullfile(dataDir, 'rgcDensityByMeridian.mat'), 'rgcDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng');
    save(fullfile(dataDir, 'rgcDisplacementMaps.mat'), 'allMaps');
end

