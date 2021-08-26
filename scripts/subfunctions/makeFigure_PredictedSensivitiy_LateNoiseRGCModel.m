function fH = makeFigure_PredictedSensivitiy_LateNoiseRGCModel(...
        prediction_visualfield, observed_visualfield, combHVA, errorCombHVA, combVMA, errorCombVMA, ...
        expName,figurePth,saveFig)

%% Bar plot to compare predictions against behavior
condNames   = {'HM', 'UVM','LVM'};
condColor   = [63, 121, 204; 0 206 209; 228, 65, 69;150 123 182]/255;
yl          = [0 4];

fH = figure(13); clf; set(fH, 'position',[139,244,1373,543], 'color', 'w', ...
    'NumberTitle', 'off', 'Name', sprintf('Predicted sensitivity given 2-AFC experiment, 4.5 deg eccen: %s', expName));
hold all;

% Plot prediction for cones
subplot(151)
bar(1:3, prediction_visualfield.cones.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(1,:)); hold on
errorbar(1:3,prediction_visualfield.cones.sensitivity.mean_wHorz,...
             (prediction_visualfield.cones.sensitivity.mean_wHorz-prediction_visualfield.cones.sensitivity.error.upper_wHorz), ...
             (prediction_visualfield.cones.sensitivity.error.lower_wHorz-prediction_visualfield.cones.sensitivity.mean_wHorz),...
             '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'linear');
set(gca, 'YLim', [0 1000], 'YScale', 'linear');
box off; ylabel('Contrast sensitivity'); title('Model Prediction Cones');

% Plot prediction for current
subplot(152)
bar(1:3, prediction_visualfield.current.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,prediction_visualfield.current.sensitivity.mean_wHorz,...
    (prediction_visualfield.current.sensitivity.mean_wHorz - prediction_visualfield.current.sensitivity.error.upper_wHorz), ...
    (prediction_visualfield.current.sensitivity.error.lower_wHorz - prediction_visualfield.current.sensitivity.mean_wHorz),...
    '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',[0 200], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'linear');
box off; ylabel('Contrast sensitivity'); title('Model Prediction Current');


% Plot prediction for mRGCs
subplot(153)
bar(1:3, prediction_visualfield.rgc.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(3,:)); hold on
errorbar(1:3,prediction_visualfield.rgc.sensitivity.mean_wHorz,...
    (prediction_visualfield.rgc.sensitivity.mean_wHorz-prediction_visualfield.rgc.sensitivity.error.upper_wHorz), ...
    (prediction_visualfield.rgc.sensitivity.error.lower_wHorz-prediction_visualfield.rgc.sensitivity.mean_wHorz), ...
    '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',[0 200], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'linear');
box off; ylabel('Contrast sensitivity'); title('Model Prediction mRGCs');

% Plot prediction for behavior from Himmelberg, Winawer, Carrasco 2020 JoV
subplot(154)
bar(1:3, observed_visualfield.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(4,:)); hold on
errorbar(1:3,observed_visualfield.sensitivity.mean_wHorz,observed_visualfield.sensitivity.error_wHorz, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',[0 200], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'linear');
box off; ylabel('Contrast sensitivity'); title('Behavior');

% Plot Asymmetries in percent
subplot(155); hold on; cla
x_bar = [0.5, 1, 1.5 2; 3 3.5 4 4.5];
for ii = 1:4
    bar(x_bar(:,ii), [combHVA(ii); combVMA(ii)], 0.2, 'EdgeColor','none','facecolor',condColor(ii,:)); hold on
end
errorbar(x_bar(1,:), combHVA, combHVA-errorCombHVA(1,:), errorCombHVA(2,:)-combHVA, '.','color', 'k', 'LineWidth',2);
errorbar(x_bar(2,:), combVMA, combVMA-errorCombVMA(1,:), errorCombVMA(2,:)-combVMA, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0,5],'Ylim',[-15 50], 'TickDir', 'out', 'XTick', [1.25, 4], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Asymmetry');

% Save figure if requested
if saveFig
    hgexport(fH, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.eps'))
    savefig(fH, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.fig'))
    print(fH, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.png'), '-dpng')
end
