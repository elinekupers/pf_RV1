function plotConeDensityIsetbio(angDeg, eccDeg, sourceDataset)
%% plotConeDensityIsetbio(angDeg, eccMM, sourceDataset)
% 
% function to plot cone density using ISETBIO toolbox
% 
% INPUTS:
%   angDeg          :   vector of polar angles to plot (in degrees)
%   eccDeg          :   vector of eccentricities to plot (in degrees)
%   [sourceDataset] :   string defining source of cone density data, choose
%                       from: 'Curcio1990' (default), 'Song2011Young' or 
%                             'Song2011Old'



% Check inputs
if ~exist('sourceDataset','var') || isempty(sourceDataset)
    sourceDataset = 'Curcio1990';
end

% Transform degrees to mm, because coneDensityReadData only takes mm as
% inputs for eccentrity values 
deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;
eccMM   = eccDeg.*deg2mm; % mm


% Allocate space
coneDensityMM2  = NaN(length(angDeg), length(eccMM));
coneDensityDeg2 = coneDensityMM2;

% Get cone density
for ii = 1:length(angDeg)
    for jj = 1:length(eccMM)
        coneDensityMM2(ii,jj) = coneDensityReadData('coneDensitySource', sourceDataset, ...
            'eccentricity', eccMM(jj), ...
            'angle', angDeg(ii), ...
            'whichEye', 'left', ...
            'angleUnits', 'deg', ...
            'eccentriticyUnits', 'mm');
    end
    
    coneDensityDeg2(ii,:) = coneDensityMM2(ii,:)./(mm2deg.^2);

end


[A, T] = meshgrid(angDeg, eccDeg);
[X, Y] = pol2cart(A, T);

fH1 = figure(); clf; set(gcf, 'Color', 'w', 'Position', [686, 345, 1223, 1000])
contour(X,Y,log10(coneDensityDeg2'))

% Add labels, make plot pretty
colormap(hsv)
xlabel('Position (deg)')
ylabel('Position (deg)')
grid on;

title('Cone density (Curcio data from ISETBIO left eye)');
c = colorbar; c.TickDirection = 'out'; ylabel(c, 'log_{10} Cones / deg^2 ');
set(gca, 'FontSize', 14', 'TickDir', 'out'); axis square;



return