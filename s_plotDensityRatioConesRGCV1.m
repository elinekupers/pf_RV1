% s_plotDensityRatioConesRGCV1


% Unit Converters
deg2m   = 0.3 * 0.001;
m2deg   = 1/deg2m;

deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Define eccentricity range
eccMM = logspace(-5,-1,20);
eccDeg = eccMM./mm2deg;

%% Set range
angDeg = (0:.02:1) * 2 * pi;  % Angle around 360 deg

%% Fill up a matrix with cone density
coneDensityMM2 = zeros(length(eccMM), length(angDeg));
coneDensityDeg2 = coneDensityMM2;
for jj = 1:length(angDeg)
    for ii = 1:length(eccDeg)
        coneDensityMM2(ii, jj) = coneDensityReadData('coneDensitySource', 'Curcio1990', ...
            'eccentricity', eccDeg(ii), ...
            'angle', angDeg(jj), ...
            'whichEye', 'left', ...
            'eccentriticyUnits', 'deg');
    end
    
    coneDensityDeg2(:,jj) = coneDensityMM2(:,jj)./(mm2deg.^2);

end


%% But the surface plot looks snazzy
[A, T] = meshgrid(angDeg, eccDeg);
[X, Y] = pol2cart(A, T);

fH1 = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [686, 345, 1223, 1000])
contour(X,Y,log10(coneDensityDeg2))

colormap(hsv)
xlabel('Position (deg)')
ylabel('Position (deg)')
grid on;

title('Cone density (left eye)');
c = colorbar; c.TickDirection = 'out'; ylabel(c, 'log_{10} Cones / deg^2 ');
set(gca, 'FontSize', 14', 'TickDir', 'out'); axis square;

