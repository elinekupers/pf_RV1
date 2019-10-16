function out = stack_interp(stack, interp_scalar)

%stack_interp takes a stack of multi-resolution images and outputs another stack of multi-resolution images of the same
%dimensions. The input stack is assumed to be created by the function optics_gaussian_stack, where the kth level of the
%stack was created by convolving some image I with a Gaussian whose standard deviation is (2^(k-1))/120 degrees. The kth
%level of the output stack will estimate the result of convolving I with a Gaussian whose standard deviation is
%(log2(interp_scalar*2^k))/120 degrees.

%REQUIRED INPUTS:
%stack = multi-resoluton stack created by optics_gaussian_stack
%interp_scalar = a scalar that changes the standard deviation of the Gaussian kernel associated with the kth level from
%                (2^(k-1))/120 degrees for the input stack to (log2(interp_scalar*2^k))/120 degrees for the output
%                stack.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE INTERPOLATED STACK

%Find the floor, ceiling and fractional interpolation levels
interp_levels = interp_scalar.*2.^(0:size(stack,3)-1);
levFl = floor(log2(interp_levels))+1;
levCl = ceil(log2(interp_levels))+1;
levFrac = log2(interp_levels)+1-levFl; 

%Find the maximum interpolation level possible
max_level = size(stack,3);
levFl(levFl>max_level) = max_level;
levCl(levCl>max_level) = max_level;

%Create stack of interpolated images
interp_stack = zeros(size(stack));
for ss = 1:size(stack,3)
    interp_stack(:,:,ss) = stack(:,:,levFl(ss)) + levFrac(ss)*(stack(:,:,levCl(ss))-stack(:,:,levFl(ss)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out = interp_stack;

