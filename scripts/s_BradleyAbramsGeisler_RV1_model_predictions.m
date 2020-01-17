%% s_BradleyAbramsGeisler_RV1_model_predictions.m

%% Define target locations (we assume uniform gray background)
% Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccentricities   = 1:9; % deg

locLabels =  {'EAST', 'NORTHEAST', 'NORTH', 'NORTHWEST', 'WEST', 'SOUTHWEST',  'SOUTH', 'SOUTHEAST'};
theta   = 0:pi/4:2*pi; % every 45 degrees, in radians

backgroundType = 'uniform'; % choose from '1f_default', '1f_recomputed' or 'uniform'
pixperdeg = 40; % pixels per 1 degree

% Plotting params
figureDir = fullfile(pfRV1rootPath, 'figures');
cmap = hsv(length(eccentricities));

for ii = 1:length(eccentricities)
    
    eccen   = eccentricities(ii); % deg
    rho     = ones(size(theta))*eccen; % deg

    [tx, ty] = pol2cart(theta, rho);

    %% Run model
    out = retina_V1_model_PF_wrapper(tx,ty, backgroundType, pixperdeg);

    %% Convert thresholds to sensitivity
    sensitivity(:,ii) = 1./out.threshold;
    
    %% Debug: plot locations of stimuli
    % figure(1); clf; set(gcf,'Color', 'w');
    % p = polarplot(theta, rho, 'ko-','LineWidth', 3);
    % title('Stimulus locations');
    % set(gca, 'FontSize', 14)
end


%% Visualize results

figure; clf; set(gcf,'Color', 'w');
for jj = 1:length(eccentricities)
   
    eccen = eccentricities(jj);
    sens  = sensitivity(:,jj);
    
    polarplot(theta, sens, 'o-', 'Color', cmap(jj,:), 'LineWidth', 3); hold on;
    
    locLabels{jj} = sprintf('Eccen %1.0f deg', eccen);
end
title('Predicted performance (contrast sensitivity)');
legend(locLabels); legend boxoff

maxr = round(max(sensitivity(:)))+1;
allRTicks = [0:10:maxr];
rlim([0 maxr])
rticks(allRTicks)
rticklabels(sprintfc('%1.0f dB',allRTicks))
    set(gca, 'FontSize', 14)

    fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
    fprintf('Predicted Horizontal-Vertical Asymmetry (sensitivity):\t %1.0f%%\n', hva(sens(1:8)))
    fprintf('Predicted Vertical-Meridian Asymmetry (sensitivity):  \t %1.0f%%\n', vma(sens(1:8)))

% Save matlab fig and pdf
figName = sprintf('PolarPlot_Sensitivity_Bradley_et_al_2014_eccen%d-%ddeg_%s_%dppd', eccentricities(1), eccentricities(end), backgroundType, pixperdeg);
savefig(fullfile(figureDir, figName))
print(fullfile(figureDir, figName), '-depsc')
print(fullfile(figureDir, figName), '-dpng')

%% Plot VMA and HVA
titleStr = sprintf('Asymmetries of sensitivity (Visual field) - Retina V1 Model %s %dppd', backgroundType, pixperdeg);
visualFieldFlag = true;
saveFigs = true;
fH1 = plotHVAandVMA(sensitivity(1:8,:), eccentricities, visualFieldFlag, titleStr, figureDir, saveFigs);


% % Print out RGC spacing
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
