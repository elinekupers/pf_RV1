%% s_BradleyAbramsGeisler_RV1_model_predictions.m

%% Define target locations (we assume uniform gray background)
% Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccen   = 9; % deg
theta   = 0:pi/4:2*pi; % every 45 degrees, in radians
rho     = ones(size(theta))*eccen; % deg

[tx, ty] = pol2cart(theta, rho);

locLabels =  {'EAST', 'NORTHEAST', 'NORTH', 'NORTHWEST', 'WEST', 'SOUTHWEST',  'SOUTH', 'SOUTHEAST'};


%% Debug: plot locations of stimuli
% figure(1); clf; set(gcf,'Color', 'w');
% p = polarplot(theta, rho, 'ko-','LineWidth', 3);
% title('Stimulus locations');
% set(gca, 'FontSize', 14)

%% Run model
out = retina_V1_model_PF_wrapper(tx,ty);

%% Convert thresholds to sensitivity
sensitivity = 1./out.threshold;

%% Visualize results

figure; 
clf; set(gcf,'Color', 'w');
%hold on
p = polarplot(theta, sensitivity, 'ko--', 'LineWidth', 3);
title(sprintf('Predicted performance (contrast sensitivity) at %d deg eccen',eccen));
set(gca, 'FontSize', 14)

hva = @(x) 100*(mean(x([1 5])) - mean(x([3 7]))) ./ mean(x([1 3 5 7]));
vma = @(x) 100*(x(7)-x(3)) / mean(x([3 7]));


fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
fprintf('Predicted Horizontal-Vertical Asymmetry (sensitivity):\t %1.0f%%\n', hva(sensitivity))
fprintf('Predicted Vertical-Meridian Asymmetry (sensitivity):  \t %1.0f%%\n', vma(sensitivity))


% Print out RGC spacing
t0 = spacing_fn(tx,ty); % RGC spacing at target center location (in deg) 

for ii = 1:length(tx)-1    
    fprintf('RGC spacing at %s: %1.3f deg\n', locLabels{ii}, t0(ii))
end    

%% RGC spacing HVA, VMA
rgcDensity = (1./t0).^2;

fprintf('Horizontal-Vertical Asymmetry (RGC density):\t %1.0f%%\n', hva(rgcDensity))
fprintf('Vertical-Meridian Asymmetry (RGC density):  \t %1.0f%%\n', vma(rgcDensity))

