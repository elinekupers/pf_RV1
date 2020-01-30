%% s_BradleyAbramsGeisler_RV1_model_predictions.m
%
% Script to run a modified version of the retina-V1 (RV1) model first
% published by Bradley, Abrams and Geisler in 2014 (JoV).
%
% Model predicts the threshold for a 2IFC detection task of a Gabor on a
% uniform or 1/f noise background. 
% We modified the model to predict thresholds for Gabors at different
% locations and to recompute the 1/f background so that wil have the same
% size as the uniform background and is not limited to the example noise 
% image.


%% PARAMETERS TO CHANGE:

% Define target locations (we assume uniform gray background)
% Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccentricities   = 1:9; % deg
locLabels        = {'EAST', 'NORTHEAST', 'NORTH', 'NORTHWEST', 'WEST', 'SOUTHWEST',  'SOUTH', 'SOUTHEAST'};
theta            = 0:pi/4:2*pi; % every 45 degrees, in radians

backgroundType   = '1f_recomputed'; % choose from '1f_default', '1f_recomputed' or 'uniform'
pixperdeg        = 40;              % pixels per 1 degree (lower numbers ~= higher SF)
nIter            = 100;             % number of simulation iterations 
verbose          = false;           % print masking coefficients or not

% Plotting params
figureDir       = fullfile(pfRV1rootPath, 'figures');
cmap            = hsv(length(eccentricities));

%% RUN BRADLEY MODEL

% Preallocate space
sensitivity =  NaN(length(theta), length(eccentricities), nIter);
T_pool      = sensitivity; % target pool response
NB          = sensitivity; % narrowband masking response
BB          = sensitivity; % broadband masking response

for ii = 1:length(eccentricities)
    
    eccen   = eccentricities(ii);       % deg
    rho     = ones(size(theta))*eccen;  % deg
    
    [tx, ty] = pol2cart(theta, rho);    % deg 
    
    for iter = 1:nIter
        
        %% Run model
        out = retina_V1_model_PF_wrapper(tx,ty, backgroundType, pixperdeg, verbose);
        
        %% Convert thresholds to sensitivity
        sensitivity(:,ii, iter) = 1./out.threshold;
        
        %% Same pooling response for target, and the narrow band/ broadband/ baseline pooling resposnes
        T_pool(:,ii, iter) = out.T_pool;
        NB(:,ii, iter) = out.NB;
        BB(:,ii, iter) = out.BB;
        
        %% Debug: plot locations of stimuli
        % figure(1); clf; set(gcf,'Color', 'w');
        % p = polarplot(theta, rho, 'ko-','LineWidth', 3);
        % title('Stimulus locations');
        % set(gca, 'FontSize', 14)
    end
    
end

%% COMPUTE MEDIAN AND STD ACROSS ITERATIONS

% Compute HVA and VMA for each eccentricity and each iteration
for iter = 1:nIter
    for ecc = 1:length(eccentricities)
        hva_tmp(iter,ecc) = hva(sensitivity(1:8,ecc,iter));
        vma_tmp(iter,ecc) = vma(sensitivity(1:8,ecc,iter));
    end
end

% Take mean and std across iterations
mdSensitivityHVA = mean(hva_tmp,1);
mdSensitivityVMA = mean(vma_tmp,1);

sdSensitivityHVA = std(hva_tmp);
sdSensitivityVMA = std(vma_tmp);

% Concatenate std errors for plotting
stdError =  vertcat(sdSensitivityHVA,sdSensitivityVMA);

%% VISUALIZE RESULTS

% Make polar plots
figure; clf; set(gcf,'Color', 'w');
for jj = 1:length(eccentricities)
    
    eccen = eccentricities(jj);
    sens  = mdSensitivity(:,jj);
    
    polarplot(theta, sens, 'o-', 'Color', cmap(jj,:), 'LineWidth', 3); hold on;
    
    locLabels{jj} = sprintf('Eccen %1.0f deg', eccen);
    
    fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
    fprintf('Predicted Horizontal-Vertical Asymmetry (sensitivity):\t %1.0f%%\n', hva(sens(1:8)))
    fprintf('Predicted Vertical-Meridian Asymmetry (sensitivity):  \t %1.0f%%\n', vma(sens(1:8)))
 
end

% Make pretty axes and labels
title('Predicted performance (contrast sensitivity)');
legend(locLabels); legend boxoff

maxr = round(max(mdSensitivity(:)))+1;
allRTicks = [0:5:maxr];
rlim([0 maxr])
rticks(allRTicks)
rticklabels(sprintfc('%1.0f dB',allRTicks))
set(gca, 'FontSize', 14)

% Save matlab fig and pdf
figName = sprintf('PolarPlot_Sensitivity_Bradley_et_al_2014_eccen%d-%ddeg_%s_%dppd_iter%d', eccentricities(1), eccentricities(end), backgroundType, pixperdeg, nIter);
savefig(fullfile(figureDir, figName))
print(fullfile(figureDir, figName), '-depsc')
print(fullfile(figureDir, figName), '-dpng')


%% PLOT VMA and HVA vs ECCENTRICITY
titleStr        = sprintf('Asymmetries of sensitivity (Visual field) - Retina V1 Model %s %dppd', backgroundType, pixperdeg);
visualFieldFlag = true; % locations are in visual field coordinates, not retinal coordinates
saveFigs        = true;

fH1 = plotHVAandVMA(mdSensitivity(1:8,:), stdError, eccentricities, visualFieldFlag, titleStr, figureDir, saveFigs);


%% Plot out RGC spacing
% t0 = spacing_fn(tx,ty); % RGC spacing at target center location (in deg)
%
% for ii = 1:length(tx)-1
%     fprintf('RGC spacing at %s: %1.3f deg\n', locLabels{ii}, t0(ii))
% end
%
% % RGC spacing HVA, VMA
% rgcDensity = (1./t0).^2;
%
% fprintf('Horizontal-Vertical Asymmetry (RGC density):\t %1.0f%%\n', hva(rgcDensity))
% fprintf('Vertical-Meridian Asymmetry (RGC density):  \t %1.0f%%\n', vma(rgcDensity))



return


%% RGC spacing HVA, VMA vs eccen
clear spacingVisualField rgcDensityVisualField
for ii = 1:length(eccentricities)
    
    eccen   = eccentricities(ii); % deg
    rho     = ones(size(theta))*eccen; % deg
    
    [tx, ty] = pol2cart(theta, rho);
    
    % Print out RGC spacing
    spacingVisualField(ii,:) = spacing_fn(tx,ty); % RGC spacing at target center location (1/sqrt(density)
    % relative to the fovea (in deg) in
    % visual field coords
    
    rgcDensityVisualField(ii,:) = (1./spacingVisualField(ii,:)).^2; % visual field coords
    
end

% Resave cardinals into retinal coords:
spacingRetina = NaN(4,length(eccentricities)); % nasal, superior, temporal, inferior RETINA)
spacingRetina(1,:) = spacingVisualField(:,rad2deg(theta) == 0);   % 1. nasal
spacingRetina(2,:) = spacingVisualField(:,rad2deg(theta) == 270); % 2. superior (flip inferior to superior)
spacingRetina(3,:) = spacingVisualField(:,rad2deg(theta) == 180); % 3. temporal
spacingRetina(4,:) = spacingVisualField(:,rad2deg(theta) == 90);  % 4. inferior (flip superior to inferior)

rgcDensityRetina = (1./spacingRetina).^2;

titleStr = 'RGC RF spacing Drasdo 2007 in retinal coords - Bradley et al 2014 - Retina V1 Model';
plotMeridiansVsEccen(spacingRetina, eccentricities, titleStr, [], true);

titleStr = 'RGC RF density Drasdo 2007 in retinal coords - Bradley et al 2014 - Retina V1 Model';
plotMeridiansVsEccen(rgcDensityRetina, eccentricities, titleStr, [], true);

titleStr = 'HVA VMA RGC RF spacing Drasdo 2007 in retinal coords - Bradley et al 2014 - Retina V1 Model';
plotHVAandVMA(spacingRetina, eccentricities, titleStr, true);

titleStr = 'HVA VMA RGC density Drasdo 2007 in retinal coords - Bradley et al 2014 - Retina V1 Model';
plotHVAandVMA(rgcDensityRetina, eccentricities, titleStr, true);
