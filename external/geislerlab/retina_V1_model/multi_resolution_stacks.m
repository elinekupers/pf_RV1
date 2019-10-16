function out = multi_resolution_stacks(target, background, Parameters)

%multi_resolution_stacks creates all the multi-resolution stacks used by retina_V1_model. The standard deviation of the
%Gaussian kernel used for the first level of each multi-resolution stack is 1/120 degrees. Pixels per degree is fixed at
%120.

%REQUIRED INPUTS:
%target = grayscale target image, range: [-127,127]. Note that the target is defined as a contrast image that can be 
%         added to the background, which has a range of [0,255]. When added to the background, [-127,127] maps to 
%         [1,255]. Thus, zero maps to 128.
%background = grayscale background image, range: [0,255].

%OPTIONAL INPUTS:
%Parameters = a structure whose fields specify the parameters of the RV1 model. To use default parameter values, use [].
%             For example: X = multi_resolution_stacks(target,background,[])

%FIELDS OF 'Parameters':
%ks = scalar that specifies the size of the surround region of a ganglion cell receptive field. Range: ks >= 1
%wc = weight on center in a Difference of Gaussians model of ganglion cell receptive fields. Range: [0,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFAULT PARAMETERS

%Default values for 'Parameters'
if isempty(Parameters)
    Parameters.ks = 10.1; %Range: ks > kc
    Parameters.wc = 0.53; %Range: [0,1]
end
P = Parameters; %rename 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE STACKS

%Make target a square image 
target = center_matrix(target,zeros(max(size(target,1),size(target,2)))); %dimensions are for smallest possible square

%Target and background Optics and Gaussian stacks
T = optics_gaussian_stack(target, 1, []); %size of target stack determined by optics_gaussian_stack
B = optics_gaussian_stack(background, 1, size(T,3)); %use size for target stack to set background stack size

%Target envelope stack
target_envelope = envelope(target);
E = optics_gaussian_stack(target_envelope, 0, size(T,3)); %do not apply optics, just convolve with Gaussians
for ee = 1:size(T,3); 
    E(:,:,ee) = E(:,:,ee)./max(max(E(:,:,ee))); %normalize envelope
end

%Difference of Gaussians stacks
T_sur = stack_interp(T, P.ks); B_sur = stack_interp(B, P.ks); %surround stacks for target and background
DT = P.wc.*T - (1-P.wc).*T_sur; %target Difference of Gaussians stack
DB = P.wc.*B - (1-P.wc).*B_sur; %background Difference of Gaussians stack
clear T_sur B_sur; %save memory

%Narrowband filter stack
NF = bandpass_filter_stacks(T); %size of NF will be [2*sx,2*sy,sz], where [sx,sy,sz] are the dimensions of T. Thus, the 
                                %size of each image in NF will be the size of a background patch. Each background patch
                                %is defined to have dimensions [2*sx,2*sy]. This convention is used because the target
                                %envelope is often larger than the input target image.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out.target_stack = T;
out.background_stack = B;
out.target_envelope_stack = E;
out.DoG_target_stack = DT;
out.DoG_background_stack = DB;
out.narrowband_filter_stack = NF;
