function stimTemplate =  getStimTemplateForSVMClassification(...
    baseFolder, subFolder, expName, cone2RGCRatio, c, eccen, selectTimePoints)
% Function is called by wrapper function linearRGCModel_ClassifyStimTemplate 
% (or absorptionsClassifierAccuracyStimTemplateWrapper) and defines stimulus 
% templates in cone or RGC space. Templates are defined from noiseless, 
% mean-subtracted absorptions for CW and CCW stimuli in quadriture phases.
% 
% INPUTS
%
% OUTPUTS
%   stimTemplate     : struct with noiseless RGC response (rows x cols) as 
%                       a template

% if ratio is empty, we assume that we deal with cone absorptions, not rgc responses
if isempty(cone2RGCRatio) || ~exist('cone2RGCRatio', 'var') 
    if strcmp(expName, 'conedensity')
        templateDir = fullfile(baseFolder, 'data', 'conedensitytemplate', 'absorptions', subFolder);
    else
        templateDir = fullfile(baseFolder, 'data', 'template', 'absorptions', subFolder);
    end
    
    templateFile = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen%1.2f_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat', ...
        c, eccen);
    
    templateNoiseless = load(fullfile(templateDir, templateFile));
    templateResponses = squeeze(templateNoiseless.absorptions(1,:,:,selectTimePoints,:));
     
else % otherwise we assume we deal with rgc responses
    
    if strcmp(expName, 'conedensity')
        templateDir = fullfile(baseFolder, 'data', 'conedensitytemplate', 'rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',cone2RGCRatio));
    else
        templateDir = fullfile(baseFolder, 'data', 'template','rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',cone2RGCRatio));
    end
    
    templateFile = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%1.2f_%s.mat', cone2RGCRatio, c, eccen, 'absorptions');
    
    templateNoiseless = load(fullfile(templateDir, templateFile));
    templateResponses = squeeze(templateNoiseless.rgcResponse(1,:,:,selectTimePoints,:));
end

    % subtract out the mean and average across time
    templateResponses = templateResponses - mean(templateResponses(:));
    templateResponses = mean(templateResponses,3);
    
    % Get rgc responses to CW and CCW, 0 and 90 phase shifted Gabors
    stimTemplate.ccw_ph1 = templateResponses(:,:,1);
    stimTemplate.ccw_ph2 = templateResponses(:,:,2);
    stimTemplate.cw_ph1  = templateResponses(:,:,3);
    stimTemplate.cw_ph2  = templateResponses(:,:,4);

return



