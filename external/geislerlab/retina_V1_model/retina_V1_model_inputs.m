function out = retina_V1_model_inputs()

%retina_V1_model_inputs calls retina_V1_model with a set of example inputs. The required and optional inputs for
%retina_V1_model are listed below.

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
target = imread('GaborPatch12.tif'); %range at this stage is [1 255]
target = cast(target,'double'); %cast to doubles
target = target - 128; %make range to [-127 127]

%A background image whose range is [0 255] at 120 pixels per degree (comment out for desired background)
background = 128*ones(4000); %uniform background
%load background_1f; background = background_1f; %1/f noise background 
%load background_natural; background = background_natural; %natural scene background 

%Fixation coordinate vectors in degrees assuming (0,0) is at center of background.
fx = [0,0,0,0];
fy = [0,0,0,0];

%Coordinate vectors of target in degrees assuming (0,0) is at the center of the background. 
tx = [0,1,2,3];
ty = [0,0,0,0];

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
Stacks = multi_resolution_stacks(target, background, Parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT (use one of the following)

%Using only required inputs:
%out = retina_V1_model(target, background, tx, ty, fx, fy, [], []);

%Using required and optional inputs:
out = retina_V1_model(target, background, tx, ty, fx, fy, Parameters, Stacks);


