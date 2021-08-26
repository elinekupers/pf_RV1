function [meshFit, gof] = ogSurfaceFit(X, Y, Z, whichfit)

% Fit model to data.
[xData, yData, zData] = prepareSurfaceData( X, Y, Z );

% Set up fittype and options.
switch whichfit
    case 'lowess'
        ft = fittype( 'lowess' );
        opts = fitoptions( 'Method', 'LowessFit' );
        opts.Normalize = 'on';
        opts.Robust = 'Bisquare';
        opts.Span = 0.2; %0.15;
        
        % Fit model to data.
        [meshFit, gof] = fit( [xData, yData], zData, ft, opts );

        
    case 'interpolant'        
        % Set up fittype and options.
        ft = 'biharmonicinterp';   
        % Fit model to data.
        [meshFit, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

end

end