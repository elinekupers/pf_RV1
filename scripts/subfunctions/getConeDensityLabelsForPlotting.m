function [labels, allDensity] = getConeDensityLabelsForPlotting(expParams)
    % Get parameters to compute cone density levels
    whichEye        = 'left';
    cparams.cmFOV   =  2; % degrees

    % Convertion deg to m
    deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter

    % Predefine density vector
    allDensity = nan(length(expParams.eccentricities),1);
    labels = cell(length(expParams.eccentricities),1);

    for ec = expParams.eccentricities
        % Specify retinal location where stimulus is presented
        cparams.eccentricity = ec;                     % Visual angle of stimulus center, in deg
        cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
        x = x * deg2m;  y = y * deg2m;
        cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

        % Set the field of view (degrees)
        cMosaic.setSizeToFOV(cparams.cmFOV);
        allDensity(ec==expParams.eccentricities,:) = eccen2density(cMosaic, 'deg');

        labels{ec==expParams.eccentricities} = sprintf('%2.2f x10^3 cells/deg2', allDensity(ec==expParams.eccentricities)/10.^3);
    end
end