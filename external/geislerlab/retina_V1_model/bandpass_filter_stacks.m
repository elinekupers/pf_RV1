function out = bandpass_filter_stacks(T)

%bandpass_filter_stacks takes a stack T of multi-resolution images and outputs a stack of filtered amplitude spectra for
%each image in T. The filters have their origin at (1,1) for use in the frequency domain. The size of the output stack
%is [2*sx,2*sy,sz], where [sx,sy,sz] are the dimensions of T. Thus, the size of each image in the output will be the
%size of a background patch. Each background patch is defined to have dimensions [2*sx,2*sy]. This convention is used
%because the target envelope is often larger than the input target image.

%Default parameters
BW_logGabor = 1.5; %radial bandwidth in octaves
BW_orient = 40/360*2*pi; %angular bandwidth in radians

%Output
szT = size(T);
Z = zeros(2*szT(1), 2*szT(2), szT(3)); Z0 = Z(:,:,1);
out = Z;

for ss = 1:size(T,3)    
    SZ = center_matrix(rot90(T(:,:,ss),2),Z0); %center target image in a matrix of zeros
    F = ifftshift(bandpass_filter(SZ, BW_logGabor, BW_orient)); %create the bandpass filter for the target image
    
    %Normalize and remove NaN
    F = F./max(abs(F(:))); F(isnan(F))=0;
    out(:,:,ss) = F;
end
    




