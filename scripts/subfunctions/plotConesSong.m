function conesSong = plotConesSong()
% Plot cone density from Song data between 1-6 degrees eccentricity.
% These data were used on Noah Benson's SFN 2019 poster.

%% Cones
conesSong = load(fullfile(pfRV1rootPath, 'external', 'data', 'song_cone_density.mat'));

% Visualize cone data
colors = {'r','b','g','k'}; % for polar angles [0 90 180 270]
colorLabels = {'temporal retina', 'superior retina', 'nasal retina', 'inferior retina'};

figure; clf; set(gcf, 'Color', 'w'); hold all;
plot(conesSong.eccentricity, conesSong.temporal, 'LineWidth',3, 'Color', colors{1});
plot(conesSong.eccentricity, conesSong.superior, 'LineWidth',3,'Color', colors{2});
plot(conesSong.eccentricity, conesSong.nasal, 'LineWidth',3,'Color', colors{3});
plot(conesSong.eccentricity, conesSong.inferior, 'LineWidth',3,'Color', colors{4});
xlabel('Eccentricity (deg)')
ylabel('Cone density (counts/deg^2)')
title('Cone density from Song et al. 2011')
set(gca, 'TickDir', 'out', 'FontSize', 20)
legend(colorLabels, 'FontSize', 20); legend boxoff;

% Assuming the delta angle for each meridian is 15 deg.
eccen = (conesSong.eccentricity>=1) & (conesSong.eccentricity<=6); % deg

integralPA15.temporal = trapz(conesSong.temporal(eccen));
integralPA15.superior = trapz(conesSong.superior(eccen));
integralPA15.nasal    = trapz(conesSong.nasal(eccen));
integralPA15.inferior = trapz(conesSong.inferior(eccen));

dataIn = [integralPA15.temporal, integralPA15.superior, integralPA15.nasal, integralPA15.inferior];

fprintf('Cone density, 1-6 deg eccen, 15 deg polar angle wedge')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(dataIn))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(dataIn))

end