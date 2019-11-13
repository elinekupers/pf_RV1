%% s_BradleyAbramsGeisler_RV1_model_predictions.m

%% Define target locations (we assume uniform gray background)
% Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccentricities   = 1:1:12; % deg

locLabels =  {'EAST', 'NORTHEAST', 'NORTH', 'NORTHWEST', 'WEST', 'SOUTHWEST',  'SOUTH', 'SOUTHEAST'};
theta   = 0:pi/4:2*pi; % every 45 degrees, in radians

% sensitivity = NaN(length(eccentricities),length(theta));

for ii = 1:length(eccentricities)
    
    eccen   = eccentricities(ii); % deg
    rho     = ones(size(theta))*eccen; % deg

    [tx, ty] = pol2cart(theta, rho);

    %% Debug: plot locations of stimuli
    % figure(1); clf; set(gcf,'Color', 'w');
    % p = polarplot(theta, rho, 'ko-','LineWidth', 3);
    % title('Stimulus locations');
    % set(gca, 'FontSize', 14)

    %% Run model
    out = retina_V1_model_PF_wrapper(tx,ty);

    %% Convert thresholds to sensitivity
    sensitivity(ii,:) = 1./out.threshold;

end



%% Visualize results
figureDir = fullfile(pfRV1rootPath, 'figures');
cmap = hsv(length(eccentricities));
 figure; clf; set(gcf,'Color', 'w');
for jj = length(eccentricities):-1:2
   
    eccen = eccentricities(jj);
    sens  = sensitivity(jj,:);
    
    polarplot(theta, sens, 'o-', 'Color', cmap(jj,:), 'LineWidth', 3); hold on;
    title('Predicted performance (contrast sensitivity) at 2-12 deg eccen');
    rlim([0 25])
    set(gca, 'FontSize', 14)

    fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
    fprintf('Predicted Horizontal-Vertical Asymmetry (sensitivity):\t %1.0f%%\n', hva(sens(1:8)))
    fprintf('Predicted Vertical-Meridian Asymmetry (sensitivity):  \t %1.0f%%\n', vma(sens(1:8)))
    
    Labels{length(eccentricities)-jj+1} = sprintf('Eccen %1.0f', eccen);
end
legend(Labels); legend boxoff

% Save matlab fig and pdf
figName = sprintf('sensitivity_Bradley_et_al_2014_eccen2-12deg');
savefig(fullfile(figureDir, figName))
print(fullfile(figureDir, figName), '-dpdf', '-fillpage')

titleStr = 'HVA VMA sensitivity Bradley et al 2014 - Retina V1 Model';
fH1 = plotHVAandVMA(sensitivity(:,1:8)', eccentricities, titleStr, true);



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




%% RGC spacing HVA, VMA vs eccen
clear t0 rgcDensity
for ii = 1:length(eccentricities)

    eccen   = eccentricities(ii); % deg
    rho     = ones(size(theta))*eccen; % deg

    [tx, ty] = pol2cart(theta, rho);

    % Print out RGC spacing
    t0(ii,:) = spacing_fn(tx,ty); % RGC spacing at target center location (in deg) 

    rgcDensity(ii,:) = (1./t0(ii,:)).^2;

end

titleStr = 'HVA VMA RGC spacing from Bradley et al 2014 - Retina V1 Model';
fH2 = plotHVAandVMA(t0(:,1:8)', eccentricities, titleStr, true);

titleStr = 'HVA VMA RGC density derived from spacing from Bradley et al 2014 - Retina V1 Model';
fH3 = plotHVAandVMA(rgcDensity(:,1:8)', eccentricities, titleStr, true);
