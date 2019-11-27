function fH = plotHVAandVMA(meridianData, regularSupportPosDegVisual, titleStr, figureDir, saveFigures)
% Function to plot HVA and VMA as a function of eccentricity. 
%
% INPUT:
%   meridianData        : cardinal meridian data (polar angle x eccen),
%                           where polar angles are in order [0 90 180 270]:
%                           (1) EAST (nasal/temporal retina, depending on eye)
%                           (2) NORTH (superior retina)
%                           (3) WEST (nasal/temporal retina, depending on eye)
%                           (4) SOUTH (inferior retina)
%                           
%                        or if both cardinal and off-cardinals in the order:
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
%   [titleStr]          : string with title for plot
%   [saveFigures]       : boolean to save figures or not
%
% OUTPUT:
%   fH                  : figure handle

% Check inputs
if isempty(titleStr) || ~exist('titleStr', 'var')
    titleStr = [];
end

if isempty(saveFigures) || ~exist('saveFigures', 'var')
    saveFigures  = false;
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
labels = {'HVA', 'VMA'};
yl     = [-80,80]; %[-1 1] * max(abs([HVAvsEccen, VMAvsEccen]));

fH = figure; clf; set(gcf, 'Color', 'w', 'Position', [809   762   751   576]); hold all;
plot(regularSupportPosDegVisual,HVAvsEccen, 'k', 'LineWidth',3);
plot(regularSupportPosDegVisual,VMAvsEccen, 'k:', 'LineWidth',3);

plot(regularSupportPosDegVisual, zeros(size(regularSupportPosDegVisual)), 'k')

h = findobj(gca, 'Type', 'line');
legend(h([3,2]), labels, 'Location', 'NorthWest'); 
legend boxoff;

ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
xlabel('Eccentricity (deg)')
title(titleStr)
set(gca', 'xlim', [0 max(regularSupportPosDegVisual)], ...
    'ylim', yl, 'TickDir', 'out', 'FontSize', 14)

if saveFigures
    % Make figure dir if doesnt exist
    if isempty(figureDir); figureDir = fullfile(pfRV1rootPath, 'figures'); end
    if ~exist(figureDir, 'dir'); mkdir(figureDir); end

    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
    
end


return