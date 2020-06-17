function [] = makeFigure1_densityVSEccentricity()

cd(pfRV1rootPath)

saveData    = false;
saveFigures = true;
loadDataFromServer = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'DissertationChapter', 'Figure1');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end


%% Define params

% Eccentricity range
eccDeg = 0:0.05:40; % deg

% Polar angle range
angDeg = [0 90 180 270]; % (nasal, superior, temporal and inferior)

xScale = 'log';
yScale = 'log';
yl = 10.^[-1 5];

%% ----------------------------------------
%  --------------- Get data ---------------
%  ----------------------------------------

% CONES from ISETBIO
conesCurcioIsetbio        = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
conesCurcioIsetbioNasal   = conesCurcioIsetbio(1,:);

% RGC from Watson
rgcWatsonIsetbio  = getMRGCRFWatson(eccDeg);
rgcWatsonIsetbioNasal = rgcWatsonIsetbio(1,:);

% V1 CMF from Horton & Hoyt
V1CMF = HortonHoytCMF(eccDeg);

%% ----------------------------------------
%  ------------ Get transforms ------------
%  ----------------------------------------

cones2mRGC_nasal = conesCurcioIsetbioNasal ./ rgcWatsonIsetbioNasal;
mRGC2V1_nasal = rgcWatsonIsetbioNasal ./ V1CMF;


%% ----------------------------------------
%  --------------- Visualize --------------
%  ----------------------------------------


fH = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [712    39   779   766])

% Cones
subplot(311)
plot(eccDeg, conesCurcioIsetbioNasal, 'k', 'LineWidth', 2);
xlim([0,max(eccDeg)]);
ylim(yl);
ylabel('Density (counts/deg^2)');
title('Cones on nasal retina (Curcio et al. 1990)')
set(gca, 'XScale', xScale, 'YScale', yScale, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
if strcmp(xScale, 'log')
     set(gca, 'XTick', [0.1 1 10 max(eccDeg)], 'XTickLabel', {'0.1' '1' '10' '40'})
end

% mRGC
subplot(312)
plot(eccDeg, rgcWatsonIsetbioNasal, 'k', 'LineWidth', 2);
xlim([0,max(eccDeg)]);
ylim(yl);
ylabel('Density (counts/deg^2)');
title('mRGC RF density on nasal retina (Watson 2015)')
set(gca, 'XScale', xScale, 'YScale', yScale, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
if strcmp(xScale, 'log')
     set(gca, 'XTick', [0.1 1 10 max(eccDeg)], 'XTickLabel', {'0.1' '1' '10' '40'})
end

% V1 CMF
subplot(313)
plot(eccDeg, V1CMF, 'k', 'LineWidth', 2);
xlim([0,max(eccDeg)]);
ylim(yl);
xlabel('Eccentricity (deg)');
ylabel('CMF (mm^2/deg^2)');
title('V1 CMF (Horton & Hoyt 1991)')
set(gca, 'XScale', xScale, 'YScale', yScale, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
if strcmp(xScale, 'log')
     set(gca, 'XTick', [0.1 1 10 max(eccDeg)], 'XTickLabel', {'0.1' '1' '10' '40'})
end

if saveFigures
    % Save matlab fig and pdf
    figName = 'Fig1A_DensityEccen';
    savefig(fH, fullfile(figureDir, figName))
    hgexport(fH, fullfile(figureDir, figName))
    print(fH, fullfile(figureDir, figName) , '-dpng')
end


%% Transformations cones to mRGC to V1
fH = figure(2); clf; set(gcf, 'Color', 'w', 'Position', [983   301   594   504])

yl2 = 10.^[-1 3];

subplot(211)
plot(eccDeg, cones2mRGC_nasal, 'k', 'LineWidth', 2);
xlim([0,max(eccDeg)]);
ylim(yl2);
xlabel('Eccentricity (deg)');
ylabel('Ratio');
title('Ratio nasal cone density : mRGC RF density')
set(gca, 'XScale', xScale, 'YScale', yScale, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','off','YMinorGrid','off','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15, 'XTick', [0.1 1 10 max(eccDeg)], 'XTickLabel', {'0.1' '1' '10' '40'});
box off;

subplot(212)
plot(eccDeg, mRGC2V1_nasal, 'k', 'LineWidth', 2);
xlim([0,max(eccDeg)]);
ylim(yl2);
xlabel('Eccentricity (deg)');
ylabel('Ratio (counts/mm^2)');
title('Ratio nasal mRGC RF density : V1 CMF')
set(gca, 'XScale', xScale, 'YScale', yScale, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','off','YMinorGrid','off','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15, 'XTick', [0.1 1 10 max(eccDeg)], 'XTickLabel', {'0.1' '1' '10' '40'});
box off;

if saveFigures
    % Save matlab fig and pdf
    figName = 'Fig1B_TransformsEccen';
    savefig(fH, fullfile(figureDir, figName))
    hgexport(fH, fullfile(figureDir, figName))
    print(fH, fullfile(figureDir, figName) , '-dpng')
end