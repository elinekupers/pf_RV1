function [] = plotCone2RGCRatioVSThreshold(expName, fit, xThresh, varargin)
% Function to plot cone density levels versus stimulus contrast thresholds.
%
% INPUTS:
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% fit             : struct with fit data
%
% xThresh         : vector with x units for plot
% [saveFig]       : boolean defining to save figure or not
% [figurePth]     : boolean defining directory where to save figure

%% 0. Parse input parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @ischar);
p.addRequired('fit', @isstruct);
p.addRequired('xThresh', @isvector);
p.addParameter('decisionmaker','SVM',  @ischar);
p.addParameter('varThresh',[],  @isvector);
p.addParameter('fitTypeName', 'linear', @ischar);
p.addParameter('yScale', 'log', @ischar);
p.addParameter('saveFig', false, @islogical);
p.addParameter('figurePth', fullfile(ogRootPath, 'figs'), @isdir);
p.parse(expName, fit, xThresh, varargin{:});

% Rename variables
expName       = p.Results.expName;
fit           = p.Results.fit;
decisionmaker = p.Results.decisionmaker;
xThresh       = p.Results.xThresh;
varThresh     = p.Results.varThresh;
fitTypeName   = p.Results.fitTypeName;
yScale        = p.Results.yScale;
saveFig       = p.Results.saveFig;
figurePth     = p.Results.figurePth;


%% 1. Fit a log-linear function to thresholds vs density

yThresh = cell2mat(fit.ctrthresh);

if strcmp(yScale,'log') % linear in log-log space
    yData = log10(yThresh);
else
    yData = yThresh;
end

if strcmp(fitTypeName, 'linear')
    [fitResult, err]  = polyfit(xThresh,yData,1);
    [y, delta] = polyval(fitResult,xThresh, err);
    R2 = 1 - (err.normr/norm(yData - mean(yData)))^2;
elseif strcmp(fitTypeName, 'poly2')
    [fitResult, err]  = polyfit(xThresh,yData,2);
    [y, delta] = polyval(fitResult,xThresh, err);
    R2 = 1 - (err.normr/norm(yData - mean(yData)))^2;
end


    
%% 2. Visualize plot
xrange = xThresh;
xticks = fliplr(xrange);
colormap = parula(length(xThresh));
for ii = 1:length(xticks); xticklbls{ii} = sprintf('%1.3f', xticks(ii)); end% get x axis range

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649]); hold on;

% plot(xThresh, 10.^y, 'r-', 'LineWidth', 3); hold all;
if ~isempty(varThresh)
    errorbar(xThresh, yThresh, varThresh, 'Color', 'k', 'LineStyle','none', 'LineWidth', 2);
end
for ii = 1:length(xThresh)
    scatter(xThresh(ii), yThresh(ii), 80, 'MarkerFaceColor', colormap(ii,:), 'MarkerEdgeColor','k', 'LineWidth',2);
end


box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,'XScale','linear', 'YScale', yScale)
xlabel('mRGC : Cone ratio','FontSize',20); ylabel('Contrast threshold (%)','FontSize',25)
set(gca, 'XTick',xticks,'XTickLabel',xticklbls, 'XLim', [0.1 2.2]);
if strcmp(yScale, 'log')
   if strcmp('SNR',decisionmaker)
       yrange = [0.01, 0.1, 1];
   elseif strcmp('SVM',decisionmaker)
       yrange = [0.001, 0.01];
   elseif strcmp('Ideal',decisionmaker)
       yrange = [0.0001, 0.001];
   end
   set(gca,'YLim', [min(yrange), max(yrange)], 'YTick',yrange,'YTickLabel',sprintfc('%1.3f',yrange.*100));
   set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
       'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
else
    set(gca, 'YTick',[0, 0.05, 0.1],'YTickLabel',[0 5 10], 'YLim', [0 0.01]);
end
legend off; title(sprintf('Contrast threshold vs mRGC:cone ratio - R^2: %1.2f', R2))

% Save fig if requested
if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('contrastThreshold_vs_Ratio%s_%s_%s_%s',expName, fitTypeName, yScale, decisionmaker)))
    hgexport(gcf,fullfile(figurePth,sprintf('contrastThreshold_vs_Ratio%s_%s_%s_%s',expName,fitTypeName, yScale, decisionmaker)))
    print(fullfile(figurePth,sprintf('contrastThreshold_vs_Ratio%s_%s_%s_%s',expName,fitTypeName, yScale, decisionmaker)), '-dpng')
end

end
