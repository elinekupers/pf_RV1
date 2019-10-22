%% s_BradleyAbramsGeisler_RV1_model_predictions.m

%Coordinate vectors of target in degrees assuming (0,0) is at the center of the background.
eccen = 6;
diag = cos(deg2rad(45))*eccen;
tx = [eccen, diag,  0,      -diag, -eccen, -diag,     0,  diag, eccen]; % horz VF (EAST), 45 deg (SOUTHEAST), LVF (SOUTH), 45 deg (SOUTHWEST), horz VF (WEST), 45 deg (NORTHWEST), UVF (NORTH), 45 deg (NORTHEAST).
ty = [0,     -diag, -eccen, -diag,      0,  diag, eccen,  diag, 0];

[theta, rho] = cart2pol(tx,ty);

% figure(1); clf; set(gcf,'Color', 'w');
% p = polarplot(theta, rho, 'ko-','LineWidth', 3);
% title('Stimulus locations');
% set(gca, 'FontSize', 14)
% 
out = retina_V1_model_PF(tx,ty);
sensitivity = 1./out.threshold;

figure(2); clf; set(gcf,'Color', 'w');
p = polarplot(theta, sensitivity, 'ko-', 'LineWidth', 3);
title(sprintf('Predicted performance (contrast sensitivity) at %d deg eccen',eccen));
set(gca, 'FontSize', 14)

hva = mean([out.threshold(3),out.threshold(7)]) / mean([out.threshold(1),out.threshold(5)]);
vma = out.threshold(7)/out.threshold(3);

fprintf('Bradley, Abrams, Geisler (2014) Retina-V1 model predictions:\n')
fprintf('Predicted Horizontal-Vertical Asymmetry:\t %1.2f\n', hva)
fprintf('Predicted Vertical-Meridian Asymmetry:  \t %1.2f\n', vma)


