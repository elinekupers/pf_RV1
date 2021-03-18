function stimTemplate =  getStimTemplateForSVMClassification(baseFolder, subFolder, expName, cone2RGCRatio, c, eccen, selectTimePoints)


    if strcmp(expName, 'conedensity')
        templateDir = fullfile(baseFolder, 'data', 'conedensitynonoise', 'rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',cone2RGCRatio));
    else
        templateDir = fullfile(baseFolder, 'data', 'template','rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',cone2RGCRatio));
    end
    
    templateFile = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%1.2f_%s.mat', cone2RGCRatio, c, eccen, 'absorptions');

    templateNoiseless = load(fullfile(templateDir, templateFile));
    templateRGCResponses = squeeze(templateNoiseless.rgcResponse(1,:,:,selectTimePoints,:));
    
    % subtract out the mean and average across time
    templateRGCResponses = templateRGCResponses - mean(templateRGCResponses(:));
    templateRGCResponses = mean(templateRGCResponses,3);
        
    % Get rgc responses to CW and CCW, 90 and 270 phase shifted Gabors
    stimTemplate.ccw_ph1 = templateRGCResponses(:,:,1);
    stimTemplate.ccw_ph2 = templateRGCResponses(:,:,2);
    stimTemplate.cw_ph1  = templateRGCResponses(:,:,3);
    stimTemplate.cw_ph2  = templateRGCResponses(:,:,4);
               
return
%

%
%     % Scene field of view
%     sceneFOV  = 2;   % scene field of view in degrees (diameter)
%     freqCPD   = 4;   % Gabor spatial frequency (cpd)
%     gausSDdeg = sceneFOV/4; % Gabor SD in degrees of visual angle
%
%     % Unit converters
%     deg2fov = 1/sceneFOV;
%     fov2deg = sceneFOV;
%
%     % Gabor parameters
%     params.freq      = fov2deg*freqCPD;
%     params.contrast  = c;                   % Presumably michelson, [0 1]
%     params.GaborFlag = gausSDdeg*deg2fov;   % Gaussian window
%     params.row       = rows;
%     params.col       = cols;
%     params.ph        = pi;  % Gabor phase (radians)
%     params1 = params;
%     params2 = params;
%     params1.ang      = (pi/180)* 15; % Gabor orientation (radians)
%     params2.ang      = -1.*params1.ang; % Gabor orientation (radians)
%
%     gaborCW_ph1      = imageHarmonic(params1);
%     gaborCCW_ph1     = imageHarmonic(params2);
%
%     params1.ph       = -pi; % Gabor phase (radians)
%     params2.ph       = -pi; % Gabor phase (radians)
%
%     gaborCW_ph2      = imageHarmonic(params1);
%     gaborCCW_ph2     = imageHarmonic(params2);
%
%     gaborCW_ph1_zeromean = gaborCW_ph1 - mean(gaborCW_ph1);
%     gaborCW_ph2_zeromean = gaborCW_ph2 - mean(gaborCW_ph2);
%
%     gaborCCW_ph1_zeromean = gaborCCW_ph1 - mean(gaborCCW_ph1);
%     gaborCCW_ph2_zeromean = gaborCCW_ph2 - mean(gaborCCW_ph2);
%
%     % Take sum of squares as a template
%     ss_ccw = ((gaborCCW_ph1_zeromean.^2) + (gaborCCW_ph2_zeromean.^2));
%     ss_cw  = ((gaborCW_ph1_zeromean.^2) + (gaborCW_ph2_zeromean.^2));
%
%     % Apply 2D FFT to  SVM template, get amplitudes
%     stimTemplate.CCW_amps = abs(fft2(ss_ccw));
%     stimTemplate.CW_amps  = abs(fft2(ss_cw));
%
%     % Set DC to 0
%     stimTemplate.CCW_amps(1,1) = 0;
%     stimTemplate.CW_amps(1,1) = 0;

%     % Normalize template
%     stimTemplate.CCW_amps = stimTemplate.CCW_amps./max(stimTemplate.CCW_amps(:));
%     stimTemplate.CW_amps = stimTemplate.CW_amps./max(stimTemplate.CW_amps(:));
%
%     % Set DC to 1
%     stimTemplate.CCW_amps(1,1) = 1;
%     stimTemplate.CW_amps(1,1) = 1;



