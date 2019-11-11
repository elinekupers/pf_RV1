function [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = plotConeAndRGCDensityDisplacementMap(sampleResEccen, sampleResPolAng, maxEccen, pixelsPerDegVisual)
%% plotConeDensityIsetbio(angDeg, eccMM, sourceDataset)
% 
% function to plot cone density using ISETBIO toolbox
% 
% INPUTS:
%   sampleResEccen      :   sampling resolution of eccentricities (in deg)
%   sampleResPolAng     :   sampling resolution of polar angles (in deg)
%   maxEccen            :   max eccentricity to plot (in deg)
%   pixelsPerDegVisual  :   number of pixels per deg visual angle, for
%                               plotting
% OUTPUTS:
%   allMaps             :   cell with cone density map {1},
%                                   rgc density map {2},
%                                   midget RGC RF {3}, 
%                                   ratio midget RGC RF to cone density {4}



% Check inputs
if ~isempty(sampleResEccen) || ~exist('sampleResEccen','var')
    sampleResEccen = 0.01;
end

if ~isempty(sampleResPolAng) || ~exist('sampleResPolAng','var')
    sampleResPolAng = 5;
end

if ~isempty(maxEccen) || ~exist('maxEccen','var')
    maxEccen = 40;
end

if ~isempty(pixelsPerDegVisual) || ~exist('pixelsPerDegVisual','var')
    pixelsPerDegVisual = 10;
end

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
     'meridianAngleResolutionDeg', sampleResPolAng); %, ...
%     'displacementMapPixelsPerDegVisual', pixelsPerDegVisual, ...
%     'coneDensityDataFileName', coneDensityDataFileName, ...
%     'rgcDensityDataFileName', rgcDensityDataFileName);


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


%% Define the image sample base for transform from polar coords
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

% loop over the maps
for vv = 1:length(polarMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(polarMapNameList{vv}), 'imRdim', imRdim);
    allMaps{vv}.data = mapImage;
    allMaps{vv}.comment = polarMapNameList{vv};
    fH1 = figure();
    fH1.Renderer='Painters';
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
    if median(arrayfun(@(x) mapImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
        opticDiscBorderColor = [0.25 0.25 0.25];
    else
        opticDiscBorderColor = [0.75 0.75 0.75];
        
    end
    plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
end

%% Plot Cone to midget RGC RF ratio
ratioConeMRFMap = allMaps{1}./allMaps{2};
inf_idx = isinf(ratioConeMRFMap);
ratioConeMRFMap(inf_idx) = 0;

fH2 = figure();
fH2.Renderer='Painters';
climVals = [0,ceil(max(ratioConeMRFMap(:)))];

displayRetinalImage(ratioConeMRFMap, climVals, pixelsPerDegVisual, maxEccen, ...
                    'Cone density / mRGCf density');
title(titleString,'FontSize',20);

%% Plot meridia cone data and midget RGC RF as separate lines

allAngles = (0:sampleResPolAng:(360-1));
cardinalMeridianAngles = [0 90 180 270];
cardinalMeridianLabels = {'nasal meridian on retina', ...
                          'superior meridian on retina', ...
                          'temporal meridian on retina', ...
                          'inferior meridian on retina'};
                          
for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(allAngles==cardinalMeridianAngles(ii));
end

meridianColors={'r','b','g','k'};

fH3 = figure(); set(gcf, 'Color', 'w', 'Position', [560   232   584   716])
clf;
for mm = 1:length(cardinalMeridianAngles)
    subplot(2,1,1); hold all;
    plot(regularSupportPosDegVisual,coneDensityByMeridian(meridianIdx(mm),:),meridianColors{mm}, 'LineWidth', 3);
    xlim([0,maxEccen]);
    ylim([0,2e4]);
    xlabel('Eccentricity (deg)');
    ylabel('Cone density [counts / deg retina ^2]');
    title('Cone density from rgcDisplacement map')
    set(gca, 'YScale', 'log', 'TickDir', 'out', 'FontSize', 14)
    if mm == 4
        legend(cardinalMeridianLabels, 'FontSize', 14); legend boxoff
    end

    subplot(2,1,2); hold all;
    plot(regularSupportPosDegVisual,mRFDensityByMeridian(meridianIdx(mm),:),meridianColors{mm}, 'LineWidth', 3);
    xlim([0,maxEccen]);
    ylim([0,3e4]);
    xlabel('Eccentricity (deg)');
    ylabel('mRF density [counts / deg retina ^2]');
    title('midget RGC RF from rgcDisplacement map')
    set(gca, 'YScale', 'log', 'TickDir', 'out', 'FontSize', 14)
    if mm == 4
        legend(cardinalMeridianLabels, 'FontSize', 14); legend boxoff
    end
end


for ii = 1:length(regularSupportPosDegVisual)
    coneHVA(ii) = hva(coneDensityByMeridian(meridianIdx,ii));
    coneVMA(ii) = vma(coneDensityByMeridian(meridianIdx, ii));
    mrgcrfHVA(ii) = hva(mRFDensityByMeridian(meridianIdx,ii));
    mrgcrfVMA(ii) = vma(mRFDensityByMeridian(meridianIdx, ii));
end

figure; clf; hold all;
plot(regularSupportPosDegVisual,coneHVA, 'k', 'LineWidth',3);
plot(regularSupportPosDegVisual,coneVMA, 'k:', 'LineWidth',3);
plot(regularSupportPosDegVisual,mrgcrfHVA, 'r', 'LineWidth',3);
plot(regularSupportPosDegVisual,mrgcrfVMA, 'r:', 'LineWidth',3);
plot(regularSupportPosDegVisual, zeros(size(regularSupportPosDegVisual)), 'k')
legend({'Cone HVA', 'Cone VMA', 'mRGC RF HVA','mRGC RF VMA', ''}, 'Location', 'Best');
legend boxoff;
ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
xlabel('Eccentricity (deg)')
title('Cone and mRGC RF density asymmetries as a function of eccentricity')
set(gca', 'xlim', [0 max(regularSupportPosDegVisual)], ...
    'ylim', [-80 80], 'TickDir', 'out', 'FontSize', 14')

