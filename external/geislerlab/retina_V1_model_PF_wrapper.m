function out = retina_V1_model_PF_wrapper(tx,ty, task, backgroundType, pixperdeg, verbose)
%
% Modification of retina_V1_model_inputs, a function that calls
% retina_V1_model with a set of example inputs.
%
% In this case, we change the polar angle of the Gabor path at iso-eccentric
% locations to see if the Retina-V1 model makes a different prediction in
% contrast thresholds.
%
% Wrapper INPUTS:
% [tx,ty]           :   (int or vector) target x,y center coordinates in visual field (deg)
% backgroundType    :   (string) define what background should be loaded
%                           choose from 'uniform' or '1f' for mean
%                           luminance or 1/f noise
% pixperdeg         :   (int) number of pixels per degree


%% Original function
%REQUIRED INPUTS:
%target = grayscale target image, range: [-127,127]. Note that the target is defined as a contrast image that can be
%         added to the background, which has a range of [0,255]. When added to the background, [-127,127] maps to
%         [1,255]. Thus, zero maps to 128.
%background = grayscale background image, range: [0,255].
%tx = x coordinate (in degrees) of target image center relative to center of background.
%ty = y coordinate (in degrees) of target image center relative to center of background.
%fx = x coordinate (in degrees) of fixation location relative to center of background.
%fy = y coordinate (in degrees) of fixation location relative to center of background.

%OPTIONAL INPUTS:
%Parameters = a structure whose fields specify the parameters of the RV1 model. To use default parameter values, use [].
%             For example: X = retina_V1_model(target,background,tx,ty,fx,fy,[],Stacks)
%Stacks = a structure whose fields consist of all the multi-resolution stacks used by the RV1 model. To create all
%         multi-resolution stacks, use the function: multi_resolution_stacks. retina_V1_model will create all necessary
%         multi-resolution stacks if 'Stacks' is not specified and the input is [].
%         For example: X = retina_V1_model(target,background,tx,ty,fx,fy,Parameters,[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REQUIRED INPUTS

%A target image whose range is [-127 127] at 120 pixels per degree
target{1} = imread('GaborPatch12.tif'); %range at this stage is [1 255]
target{1} = cast(target{1},'double'); %cast to doubles
target{1} = target{1} - 128; %make range to [-127 127]

%A background image whose range is [0 255] at 120 pixels per degree (comment out for desired background)
switch backgroundType
    case 'uniform'
        background = 128*ones(4000); %uniform background
        
    case '1f_default'
        load background_1f; background = background_1f; %1/f noise background
        
    case '1f_recomputed'
        tmp = noiseonf(4000, 1);
        tmp = tmp + abs(min(tmp(:)));
        background = (tmp./max(tmp(:)));
        contrastBounds = prctile(background(:),[2.5, 97.5]);
        background(background<contrastBounds(1)) = contrastBounds(1);
        background(background>contrastBounds(2)) = contrastBounds(2);
        
        background = background-min(background(:));
        background = 256.*(background./max(background(:)));
        
    case 'natural'
        load background_natural; background = background_natural; %natural scene background
end

if strcmp(task, '2AFC')   
    % rotate image this hacky way. Future computation should recompute
    % Gabor with params below
    
    CCW = imrotate(target{1}, -15, 'bilinear','crop');
    CW  = imrotate(target{1}, 15, 'bilinear','crop');
    
    target{1} = CCW;
    target{2} = CW;
    
%     stimparams.row = 64;
%     stimparams.col = 64;
%     stimparams.contrast = 1;
%     stimparams.ph = pi / 1;
%     stimparams.freq = 4;
%     stimparams.ang = pi/2 - (pi/6);
%     stimparams.GaborFlag = 0.15;
%     [img, p] = imageHarmonic(stimparams);
%     figure;
%     imagesc(img);
%     colormap(gray);
%     axis image
    
end
    

if verbose
    % Plot target on background
    figure; clf; 
    imagesc(background); hold all; colormap gray;
    [w, h] = size(background);
    wC = round(w/2);
    hC = round(h/2);
    for ii = 1:length(tx)
        plot(wC+(tx(ii).*pixperdeg),hC+(ty(ii).*pixperdeg), 'ko',  'MarkerSize',15, 'LineWidth', 3)
    end
    title('Stimulus locations');
    xlabel('Pixels'); ylabel('Pixels'); box off;
    set(gca, 'TickDir', 'out', 'FontSize', 12)
end
    
%Fixation coordinate vectors in degrees assuming (0,0) is at center of background.
fx = zeros(size(tx));
fy = zeros(size(ty));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OPTIONAL INPUTS

%Parameters
Parameters.sl = 120; %Range: sl >= 1
Parameters.kc = 1; %Range: kc >= 1
Parameters.ks = 10.1; %Range: ks > kc
Parameters.wc = 0.53; %Range: [0,1]
Parameters.kb = 24; %Range: kb >= 0
Parameters.wb = 0.925; %Range: [0,1]
Parameters.ro = 2.4; %Range: p0 > 0
Parameters.p0 = 0.0014; %0.0014 for ModelFest, 0.00045 for our yes-no experiment
Parameters.dp_crit = 1.8009; %corresponds to 0.5+0.5*(1-exp(-1)) = 81.61 percent

%Multi-resolution stacks
if strcmp(task, '2AFC')
    Stacks{1}  = multi_resolution_stacks(target{1}, background, Parameters);
    Stacks{2}  = multi_resolution_stacks(target{2}, background, Parameters);
else
    Stacks{1} = multi_resolution_stacks(target, background, Parameters);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT (use one of the following)

%Using only required inputs:
%out = retina_V1_model(target, background, tx, ty, fx, fy, [], []);

%Using required and optional inputs:
out = retina_V1_model_PF(target, background, task, tx, ty, fx, fy, pixperdeg, Parameters, Stacks, verbose);


