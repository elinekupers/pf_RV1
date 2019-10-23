%% s_BradleyAbramsGeisler_RV1_model_predictions.m

%% Define target locations (we assume uniform gray background)
% Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccen = 6; % deg
diag  = cos(deg2rad(45))*eccen; % deg
tx    = [eccen, diag,  0,      -diag, -eccen, -diag,     0,  diag, eccen]; % horz VF (EAST), 45 deg (SOUTHEAST), LVF (SOUTH), 45 deg (SOUTHWEST), horz VF (WEST), 45 deg (NORTHWEST), UVF (NORTH), 45 deg (NORTHEAST).
ty    = [0,     -diag, -eccen, -diag,      0,  diag, eccen,  diag, 0];

locLabels =  {'EAST','SOUTHEAST', 'SOUTH', 'SOUTHWEST', 'WEST', 'NORTHWEST', 'NORTH', 'NORTHEAST'};

[theta, rho] = cart2pol(tx,ty);

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

figure(2); clf; set(gcf,'Color', 'w');
p = polarplot(theta, sensitivity, 'ko-', 'LineWidth', 3);
title(sprintf('Predicted performance (contrast sensitivity) at %d deg eccen',eccen));
set(gca, 'FontSize', 14)

hva = mean([sensitivity(1),sensitivity(5)]) / mean([sensitivity(3),sensitivity(7)]);
vma = sensitivity(3)/sensitivity(7);

fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
fprintf('Predicted Horizontal-Vertical Asymmetry (sensitivity):\t %1.2f\n', hva)
fprintf('Predicted Vertical-Meridian Asymmetry (sensitivity):  \t %1.2f\n', vma)


% Print out RGC spacing
for ii = 1:length(tx)
    
    t0 = spacing_fn(tx(ii),ty(ii)); % RGC spacing at target center location (in deg) 
    fprintf('RGC spacing at %s: %1.3f deg\n', locLabels{ii}, t0)
end    
