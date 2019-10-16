function out = resmap_interp(stack, resmap)

%resmap_interp takes a stack of multi-resolution images and outputs an interpolated image based on resmap. It is assumed
%that the standard deviation of the Gaussian kernel associated with the first level of the stack is 1/120th of a degree,
%while the standard deviation of the Gaussian kernel for the kth level is 2^(k-1) times that for the kth level. Pixels
%per degree is fixed at 120.

%REQUIRED INPUTS:
%stack = multi-resoluton stack created by optics_gaussian_stack
%resmap = a foveated spatial resolution map created by create_resmap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE INTERPOLATED RESMAP

%Standard deviation for first level of stack is fixed at 1/120th of a degree (pixels per degree = 120)
std_init = 1/120;

%Find the floor, ceiling and fractional levels
levFl = floor(log2(resmap./std_init))+1;
levCl = levFl+1;
levFrac = log2(resmap./std_init)+1-levFl;

%Find the maximum interpolation level possible
sz = size(stack);
if length(sz)==3
    max_level = sz(3);
else
    max_level = 1;
end
levFl(levFl>max_level) = max_level;
levCl(levCl>max_level) = max_level;

%Create interpolated image
[X,Y] = meshgrid(1:sz(1),1:sz(2));
SubFl = sub2ind(sz,X',Y',levFl); SubCl = sub2ind(sz,X',Y',levCl);
resFl = stack(SubFl); resCl = stack(SubCl);
resFrac = (1-levFrac).*(resFl - resCl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out = resCl + resFrac;
