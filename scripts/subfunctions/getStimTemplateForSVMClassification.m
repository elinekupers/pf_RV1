function stimTemplate =  getStimTemplateForSVMClassification(baseFolder, subFolder, cone2RGCRatio, c, eccen, selectTimePoints)

templateDir = fullfile(baseFolder, 'data', 'conedensitynonoise', 'rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',cone2RGCRatio));
templateFile = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%1.2f_%s.mat', cone2RGCRatio, c, eccen, 'absorptions');

templateNoiseless = load(fullfile(templateDir, templateFile));
templateRGCResponses = templateNoiseless.rgcResponse;

% Get rgc responses to CW and CCW, 90 and 270 phase shifted Gabors
ccw_ph1 = templateRGCResponses(:,:,:,selectTimePoints,1);
cw_ph1  = templateRGCResponses(:,:,:,selectTimePoints,3);
ccw_ph2 = templateRGCResponses(:,:,:,selectTimePoints,2);
cw_ph2  = templateRGCResponses(:,:,:,selectTimePoints,4);

% Take mean across time points
ccw_ph1_mn = mean(ccw_ph1,4);
cw_ph1_mn  = mean(cw_ph1,4);
ccw_ph2_mn = mean(ccw_ph2,4);
cw_ph2_mn  = mean(cw_ph2,4);

% Take mean across trials
ccw_ph1_mn2 = squeeze(mean(ccw_ph1_mn,1));
cw_ph1_mn2  = squeeze(mean(cw_ph1_mn,1));
ccw_ph2_mn2 = squeeze(mean(ccw_ph2_mn,1));
cw_ph2_mn2  = squeeze(mean(cw_ph2_mn,1));

% Take sum of squares as a template
ss_ccw = ((ccw_ph1_mn2.^2) + (cw_ph1_mn2.^2));
ss_cw  = ((ccw_ph2_mn2.^2) + (cw_ph2_mn2.^2));

% Apply 2D FFT to  SVM template, get amplitudes
stimTemplate.CCW_amps = abs(fft2(ss_ccw));
stimTemplate.CW_amps  = abs(fft2(ss_cw));


return
% 
% if c == 0
%     stimTemplate.CCW_amps = ones(rows,cols);
%     stimTemplate.CW_amps = ones(rows,cols);
% else
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



