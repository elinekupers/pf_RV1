function out = optics_gaussian_stack(img, optics, max_level)

%optics_gaussian_stack creates a stack of multi-resolution images from the input image img. For each level of the stack,
%img is convolved with a Gaussian kernel. The standard deviation of this Gaussian kernel is 1/120th of a degree for the
%first level, and then is 2^(k-1) times that for the kth level. Pixels per degree is fixed at 120. If one wishes to also
%apply the effect of human optics to every level of the stack, set optics = 1. Otherwise, set optics = 0.

%The number of levels in the stack is determined by 'max_level'. 'max_level' can be set to any non-negative integer, or
%it can be set to []. If max_level = [], then the number of levels depends on properties of img. Specifically, if
%max_level = [], then the last level of the stack is the first level with variance = 0. Because it is rare for a stack
%to require more than 12 levels, setting max_level = [] creates a stack with at most 12 multi-resolution levels.

%REQUIRED INPUTS:
%img = grayscale image, range: [0,255] (e.g. 8-bit image, though real numbers in this range are acceptable)
%optics = controls whether or not to apply the effect of human optics to all levels of the multi-resolution stack. To
%         apply the optics to all levels, set optics = 1. Otherwise, set optics = 0. 

%OPTIONAL INPUTS:
%max_level = specify maximum number of levels for the stack. If specified, max_level must be a non-negative integer. To
%            use the default, set max_level = []. The default sets the number of levels so that the last level has a
%            variance = 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFAULT PARAMETERS

ppd = 120; %Pixels per degree
g_stds = (2.^(0:11))./ppd; %standard deviations of Gaussian kernels for different stack levels

%Maximum number of levels
if isempty(max_level)
    max_level = length(g_stds); %max default length = 12
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE STACK

%Optics of the human eye
sz = size(img); degX = sz(1)/ppd; degY = sz(2)/ppd; [Y,X] = meshgrid(1:sz(2),1:sz(1));
freq = sqrt(((X-center(sz(1))).^2)./(degX.^2)+((Y-center(sz(2))).^2)./(degY.^2)); %input frequencies
MTF = MTF_optics(freq); %MTF_optics allows an entire array of frequencies as input

%Fourier transform of input image
fft_img = fft2(ifftshift(img)); 

MRL = zeros(sz(1),sz(2),max_level); 
gg = 1;
while gg <= max_level
    
    %Convolve img with a Gaussian of specified standard deviation
    Gaus = fft_Gaus(g_stds(gg),sqrt(((X-center(sz(1))).^2)./(degX.^2)+((Y-center(sz(2))).^2)./(degY.^2)));
    
    %Apply optics if desired
    if optics == 1
        Q = fftshift(ifft2(fft_img.*ifftshift(MTF).*ifftshift(Gaus)));
    elseif optics == 0
        Q = fftshift(ifft2(fft_img.*ifftshift(Gaus)));
    end
    MRL(:,:,gg) = Q;
    
    %Truncate array if current level has variance = 0
    if var(Q(:)) == 0
        MRL = MRL(:,:,1:gg);
        gg = max_level;
    end
    
    gg = gg+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out = MRL;