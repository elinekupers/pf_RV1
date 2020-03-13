function out = retina_V1_model_PF_wrapper(tx,ty, contrast, sparams)
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
% contrast          :   (double) target contrast level in % between [0,1]
% sparams           :   (struct) with stimulus variables, for example:
%                         - backgroundType: (string) define what background should be loaded
%                                       choose from 'uniform' or '1f' for mean
%                                       luminance or 1/f noise
%                         - pixperdeg:  (int) number of pixels per degree
%                         - task:       (string) definition of
%                                       psychophysical task 


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
if ~exist('contrast', 'var') || isempty(contrast)
    contrast = 1;
end

switch sparams.stimType
    case 'Gabor_default'
        %A target image whose range is [-127 127] at 120 pixels per degree
        target{1} = imread('GaborPatch12.tif'); %range at this stage is [1 255]
        target{1} = cast(target{1},'double'); %cast to doubles
        target{1} = target{1} - 128; %make range to [-127 127]
        
        if strcmp(sparams.task, '2AFC')
            % rotate image this hacky way. Future computation should recompute
            % Gabor with params below
            CCW = imrotate(target{1}, (90-15), 'bilinear','crop');
            CW  = imrotate(target{1}, (90+15), 'bilinear','crop');
            target{1} = CCW;
            target{2} = CW;
        end
    case 'Gabor_observermodel'
        load(fullfile(sparams.stimDir, sprintf('stimulus_contrast%1.4f_pa0_eccen4.50_defocus0.00_sf4.00.mat',contrast)), 'scenes');
        stim      = scenes{1}{2}.data.photons(:,:,1);
        stimNorm  = stim./max(stim(:));
        target{1} = 256.*(stimNorm-prctile(stimNorm(:),50));
        target{1} = cast(target{1},'double'); %cast to doubles
        
        if strcmp(sparams.task, '2AFC')
            CW        = target{1};
            stim      = scenes{3}{2}.data.photons(:,:,1);
            stimNorm  = stim./max(stim(:));
            CCW       = 256.*(stimNorm-prctile(stimNorm(:),50));
            CCW       = cast(CCW,'double'); %cast to doubles
            
            target{1} = CCW;
            target{2} = CW;
        end
end

%A background image whose range is [0 255] at 120 pixels per degree (comment out for desired background)
switch sparams.backgroundType
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



if sparams.verbose
    % Plot target on background
    figure; clf; colormap gray; colorbar;
    [w, h] = size(background);
    wC = round(w/2);
    hC = round(h/2);
    tmp = background;
    for ii = 1:length(tx)
        stimrows = (wC + round(tx(ii).*sparams.pixperdeg)) : ((wC + round(tx(ii).*sparams.pixperdeg))+ size(target{1},2)-1);
        stimcols = (hC + round(ty(ii).*sparams.pixperdeg)) : ((hC + round(ty(ii).*sparams.pixperdeg))+ size(target{1},1)-1);
        tmp(stimrows,stimcols) = tmp(stimrows,stimcols)+target{1};
    end
    imagesc(tmp)
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
if strcmp(sparams.task, '2AFC')
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
out = retina_V1_model_PF(target, background, sparams.task, tx, ty, fx, fy, sparams.pixperdeg, Parameters, Stacks, sparams.verbose);


