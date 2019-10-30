% s_plotDensityRatioConesRGCV1

%% Define params
% Unit Converters
deg2m   = 0.3 * 0.001;
m2deg   = 1/deg2m;

deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Define eccentricity range
eccDeg = linspace(0,30,73); % deg
eccMM  = eccDeg.*deg2mm; % mm

% Define polar angle range
angDeg = (0:5:360);  % deg
% angDeg = deg2rad([0, 90, 180, 270, 360]);

%% -------- CONES ----------
clear coneDensityMM2 coneDensityDeg2

% Get cone density
for ii = 1:length(angDeg)
    for jj = 1:length(eccMM)
        coneDensityMM2(ii,jj) = coneDensityReadData('coneDensitySource', 'Curcio1990', ...
            'eccentricity', eccMM(jj), ...
            'angle', angDeg(ii), ...
            'whichEye', 'left', ...
            'angleUnits', 'deg', ...
            'eccentriticyUnits', 'mm');
    end
    
    coneDensityDeg2(ii,:) = coneDensityMM2(ii,:)./(mm2deg.^2);

end


%% Plot filled contour map
% nanMask = ~isnan(coneDensityDeg2); 

[X, Y] = pol2cart(deg2rad(angDeg),eccDeg);
[XX, YY] = meshgrid(angDeg, eccDeg);

fH1 = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [686, 345, 1223, 1000])
contourf(XX,YY,log10(coneDensityDeg2))

colormap(hsv)
xlabel('Position (deg)')
ylabel('Position (deg)')
grid on;

title('Cone density (left eye)');
c = colorbar; c.TickDirection = 'out'; ylabel(c, 'log_{10} Cones / deg^2 ');
set(gca, 'FontSize', 14', 'TickDir', 'out'); axis square;


%% -------- RGC ISETIO (under construction) ----------

% using isetbio, we first need to create a cone mosaic object
fov = 2; % deg
cMosaic = coneMosaic('center', [0 0]);
cMosaic.setSizeToFOV(fov);

% Create bipolar layer
clear bpL bpMosaicParams
bpL = bipolarLayer(cMosaic);

bpMosaicParams.spread  = 2;  % ???RF diameter w.r.t. input samples 
bpMosaicParams.stride  = 2;  % ???RF diameter w.r.t. input samples
bpL.mosaic{1} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);

% Create a rgc layer
clear rgcL rgcParams
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 2;

% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget');
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget',rgcParams);


% nr of cells
[x_cells, y_cells, ~] = size(rgcL.mosaic{1}.cellLocation);
nrOfCells = prod([x_cells,y_cells]);

% what's the size of the rgc layer?
patchSize  = prod(rgcL.size);

% Density within fov
rgcDensity = nrOfCells * (1/patchSize);

