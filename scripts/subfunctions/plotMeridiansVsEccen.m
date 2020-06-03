function fH = plotMeridiansVsEccen(meridianData, regularSupportPosDegVisual, ...
                                    titleStr, yl, figureDir, saveFigures)
% Function to plot meridia data as a function of eccentricity. 
%
% INPUT:
%   meridianData        : cardinal meridian data (polar angle x eccen),
%                           where polar angles are in order [0 90 180 270]:
%                           (1) EAST (nasal/temporal retina, depending on eye)
%                           (2) NORTH (superior retina)
%                           (3) WEST (nasal/temporal retina, depending on eye)
%                           (4) SOUTH (inferior retina)
%                           
%   regularSupportPosDegVisual : vector with eccentricity in deg of visual
%                                angle.
%   [titleStr]          : string with title for plot
%   [yl]                : vector defining y axis limits
%   [saveFigures]       : boolean to save figures or not
%
% OUTPUT:
%   fH                  : figure handle

%% Plot meridia data as a function of eccentricity
if isempty(titleStr) || ~exist('titleStr', 'var')
    titleStr = [];
end

if isempty(yl) || ~exist('yl', 'var')
    yl       = [0 max(meridianData(:))+0.1*max(meridianData(:))];
end

if isempty(saveFigures) || ~exist('saveFigures', 'var')
    saveFigures  = false;
end

% Get max eccentricity for xlim
maxEccen = max(regularSupportPosDegVisual);

% get labels and color for plotting
cardinalMeridianAngles = [0 90 180 270];
cardinalMeridianLabels = {'nasal meridian on retina', ...
                          'superior meridian on retina', ...
                          'temporal meridian on retina', ...
                          'inferior meridian on retina'};

meridianColors         = {'r','b','g','k'};

% Visualize
fH = figure(); clf; set(gcf, 'Color', 'w', 'Position', [560   421   584   527])

for mm = 1:length(cardinalMeridianAngles)
    hold all;
    plot(regularSupportPosDegVisual, meridianData(mm,:), meridianColors{mm}, 'LineWidth', 3);
    xlim([0,maxEccen]);
    ylim(yl);
    xlabel('Eccentricity (deg)');
    ylabel('Density [counts / deg retina ^2]');
    title(titleStr)
    set(gca, 'YScale', 'log', 'TickDir', 'out', 'FontSize', 14)
    if mm == 4
        legend(cardinalMeridianLabels, 'FontSize', 14); legend boxoff
    end

end

if saveFigures
    % Make figure dir if doesnt exist
    if isempty(figureDir); figureDir = fullfile(pfRV1rootPath, 'figures'); end
    if ~exist(figureDir, 'dir'); mkdir(figureDir); end

    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
    
end


