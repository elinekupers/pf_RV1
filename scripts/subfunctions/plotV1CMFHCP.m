function v1CMF = plotV1CMFHCP()
% Plot V1 cortical magnification factor (CMF, mm2 of cortex / deg visual field) 
% between 1 and 6 degrees eccentricity.
% These data were used on Noah Benson's SFN 2019 poster.

%% V1 CMF
v1CMF = readtable(fullfile(pfRV1rootPath, 'external', 'data', 'HV-DV-surface-areas.csv'));

eccen    = [1, 6]; % deg
eccenIdx = v1CMF.min_eccen==eccen(1); % deg
polang   = v1CMF.delta_angle==30; % deg

% TO CHECK
temporal = (v1CMF.meridian=={'horizontal'}) & (v1CMF.label=='v1d');
nasal    = (v1CMF.meridian=='horizontal') & (v1CMF.label=='v1v');
inferior = (v1CMF.meridian=='vertical') & (v1CMF.label=='v1d');
superior = (v1CMF.meridian=='vertical') & (v1CMF.label=='v1v');

eccentricity = linspace(eccen(1),eccen(2),1000);

% Visualize cone data
colors = {'r','b','g','k'}; % for polar angles [0 90 180 270]
colorLabels = {'temporal V1', 'superior V1', 'nasal V1', 'inferior V1'};

figure; clf; set(gcf, 'Color', 'w'); hold all;
plot(eccentricity, temporal, 'LineWidth',3, 'Color', colors{1});
plot(eccentricity, superior, 'LineWidth',3,'Color', colors{2});
plot(eccentricity, nasal, 'LineWidth',3,'Color', colors{3});
plot(eccentricity, inferior, 'LineWidth',3,'Color', colors{4});
xlim([0 7])
xlabel('Eccentricity (deg)')
ylabel('V1 CMF (mm2/deg^2)')
title('mRGC RF density from Watson 2014')
set(gca, 'TickDir', 'out', 'FontSize', 20, 'YScale', 'log')
legend(colorLabels, 'FontSize', 20); legend boxoff;

integralPA15.temporal = trapz(rgcWatson.temporal(eccen));
integralPA15.superior = trapz(rgcWatson.superior(eccen));
integralPA15.nasal    = trapz(rgcWatson.nasal(eccen));
integralPA15.inferior = trapz(rgcWatson.inferior(eccen));

dataIn = [integralPA15.temporal, integralPA15.superior, integralPA15.nasal, integralPA15.inferior];

fprintf('V1 CMF 1-6 deg eccen, 15 deg polar angle wedge')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(dataIn))
fprintf('Vertical-Meridian Asymmetry (North / South):\t %1.2f%%\n', vma(dataIn))

end


