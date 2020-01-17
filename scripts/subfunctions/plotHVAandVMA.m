function fH = plotHVAandVMA(meridianData, regularSupportPosDegVisual, visualFieldFlag, titleStr, figureDir, saveFigures)
% Function to plot HVA and VMA as a function of eccentricity. 
%
% INPUT:
%   meridianData         : cardinal meridian data (polar angle x eccen),
%                           where polar angles are in order [0 90 180 270]:
%                           (1) EAST (nasal/temporal retina, depending on eye)
%                           (2) NORTH (superior retina)
%                           (3) WEST (nasal/temporal retina, depending on eye)
%                           (4) SOUTH (inferior retina)
%                           
%                           or if both cardinal and off-cardinals in the order:
%                           (1) EAST  (nasal/temporal retina, depending on eye)
%                           (2) NORTHEAST 
%                           (3) NORTH (superior retina) 
%                           (4) NORTHWEST,
%                           (5) WEST  (nasal/temporal retina, depending on eye) 
%                           (6) SOUTHWEST, 
%                           (7) SOUTH (inferior retina)
%                           (8) SOUTHEAST
%   regularSupportPosDegVisual : vector with eccentricity in deg of visual
%                                angle.
%   [visualFieldFlag]   : boolean to note if coordinates are not in
%                           retinal coordinates but in visual field coordinates
%                           default is 0 (retinal coordinates)
%   [titleStr]          : string with title for plot
%   [figureDir]         : string with directory to save figures
%   [saveFigures]       : boolean to save figures or not
%
% OUTPUT:
%   fH                  : figure handle

% Check inputs
if isempty(visualFieldFlag) || ~exist('visualFieldFlag', 'var')
    visualFieldFlag = false;
end

if isempty(titleStr) || ~exist('titleStr', 'var')
    titleStr = [];
end

if isempty(figureDir) || ~exist('figureDir', 'var')
    figureDir  = fullfile(pfRV1rootPath, 'figures');
end

if isempty(saveFigures) || ~exist('saveFigures', 'var')
    saveFigures  = false;
end

if visualFieldFlag
    if size(meridianData,1)==4
        meridianDataVisualField = NaN(size(meridianData));
        
        % Flip East and West
        meridianDataVisualField(1,:) = meridianData(3,:);
        meridianDataVisualField(3,:) = meridianData(1,:); 
        
        % Flip North and South
        meridianDataVisualField(2,:) = meridianData(4,:);
        meridianDataVisualField(4,:) = meridianData(2,:); 
    
    elseif size(meridianData,1)==8
        meridianDataVisualField = NaN(size(meridianData));
        
        % Flip East and West
        meridianDataVisualField(1,:) = meridianData(5,:);
        meridianDataVisualField(5,:) = meridianData(1,:); 
        
        % Flip North and South
        meridianDataVisualField(3,:) = meridianData(7,:);
        meridianDataVisualField(7,:) = meridianData(3,:); 
        
        % Flip off-cardinals
        meridianDataVisualField(2,:) = meridianData(6,:);
        meridianDataVisualField(6,:) = meridianData(2,:);
        meridianDataVisualField(4,:) = meridianData(8,:);
        meridianDataVisualField(8,:) = meridianData(4,:); 
        
    end
    
    meridianData = meridianDataVisualField;
    yLabelHVA = 'higher performance on vertical VF <- Asymmetry (%) -> higher performance on horizontal VF';
    yLabelVMA = 'higher performance on upper VF <- Asymmetry (%) -> higher performance on lower VF';
    
else
    yLabelHVA = 'more vertical retina <- Asymmetry (%) -> more horizontal retina';
    yLabelVMA = 'more inf retina <- Asymmetry (%) -> more sup retina';   
end


% Define HVA and VMA for every eccentricity data point
for ii = 1:length(regularSupportPosDegVisual)
    
    % Remove blindspot data (if there)
    if ((meridianData(1,ii) == 0) && (regularSupportPosDegVisual(ii)>12.8) && (regularSupportPosDegVisual(ii)<18.05))
        
        HVAvsEccen(ii) = NaN;
        VMAvsEccen(ii) = NaN;
    else
        HVAvsEccen(ii) = hva(meridianData(:,ii));
        VMAvsEccen(ii) = vma(meridianData(:,ii));
    end
end




% Go plot
yl     = [-80,80]; %[-1 1] * max(abs([HVAvsEccen, VMAvsEccen]));

fH = figure; clf; set(gcf, 'Color', 'w', 'Position', [418, 269, 1905, 872]); hold all;
subplot(121);
plot(regularSupportPosDegVisual,HVAvsEccen, 'ko-', 'LineWidth',3); hold on;
plot(regularSupportPosDegVisual, zeros(size(regularSupportPosDegVisual)), 'k')
grid on; axis square; box off;
ylabel(yLabelHVA)
xlabel('Eccentricity (deg)')
title(['HVA ' titleStr])
set(gca', 'xlim', [0 max(regularSupportPosDegVisual)], ...
    'ylim', yl, 'TickDir', 'out', 'FontSize', 14)

subplot(122);
plot(regularSupportPosDegVisual,VMAvsEccen, 'ko-', 'LineWidth',3); hold on;
plot(regularSupportPosDegVisual, zeros(size(regularSupportPosDegVisual)), 'k')
grid on; axis square; box off;
ylabel(yLabelVMA)
xlabel('Eccentricity (deg)')
title(['VMA ' titleStr])
set(gca', 'xlim', [0 max(regularSupportPosDegVisual)], ...
    'ylim', yl, 'TickDir', 'out', 'FontSize', 14)

if saveFigures
    % Make figure dir if doesnt exist
    if ~exist(figureDir, 'dir'); mkdir(figureDir); end

    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH, fullfile(figureDir, sprintf('%s', figName)))
    print(fullfile(figureDir, figName), '-depsc')
    print(fullfile(figureDir, figName), '-dpng')
    
end


return