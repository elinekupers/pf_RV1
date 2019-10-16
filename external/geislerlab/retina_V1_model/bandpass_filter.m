function out = bandpass_filter(img, BW_logGabor, BW_orient)

%bandpass_filter outputs a bandpass filter in Fourier space for a img whose dimensions are nxn (square). The filter is
%created through a frequency space variant convolution with a Gaussian whose standard deviations change as a function of
%frequency and orientation.  The Gaussian is defined in polar coordinates with BW_logGabor specifying the bandwidth in
%octaves in the radial direction and BW_orient specifying the bandwidth in radians in angular direction.

%Suggested defaults: BW_logGabor = 1.5, and BW_orient = 40/360*2*pi.

%Since convolution wraps around, we enlarge img to prevent edge artifacts
k_scalar = 3;
max_size = k_scalar*size(img,1);
new_img = center_matrix(img,zeros(max_size));

%We take the Fourier transform of the img.
fft_img = fft2(ifftshift(new_img));
Z = abs(fft_img); %this is the amplitude spectrum of img
sz1 = size(Z,1); sz2 = size(Z,2);
cz = center(sz1);
[X,Y] = meshgrid(-floor(sz1/2):ceil(sz1/2)-1,-floor(sz2/2):ceil(sz2/2)-1);
Z_shift = ifftshift(Z);

%These are the points we want to interpolate at. The matrix size is the
%same as that of the img.
Tmax = pi+0.001; %we wrap around along theta.
newT = meshgrid(linspace(-Tmax,Tmax,sz1));
newT = newT';
maxR = log(cz-1);
newR = meshgrid(linspace(0,1,sz1)*(0.999*maxR));

%We find the coordinates in X-Y space of the interpolation points.
[X1,Y1] = pol2cart(newT,exp(newR));

%The values at the interpolation points.
Z1_shift = interp2(X,Y,Z_shift,X1,Y1);

%We need to know the conversions for the axes of Z1 in LG-VM space.
LG_conversion = (cz-1)/maxR; %pixels per log unit.
OR_conversion = sz1/(2*Tmax); %pixels per radian.

%Change bandwidths to standard deviations. Note: FWHM = 2*sqrt(2*log(2))*std.
OR_std = BW_orient/(2*sqrt(2*log(2)));
LG_std = BW_logGabor/(2*sqrt(2*log(2)));

%We define the Gaussian that we convolve with.
sigma_x = OR_std*OR_conversion;
sigma_y = LG_std*LG_conversion;
G = gaussian_fn(sigma_y,sigma_x,sz1,sz2);

%Apply proper scalar using change of variables theorem
Z1_shift1 = Z1_shift.*exp(newR);

%Stretch so that edges aren't artificially altered during convolution
Z1L = Z1_shift1(:,1); Z1R = Z1_shift1(:,end);
Z1_shift2(:,1:sz2) = repmat(Z1L,1,sz2);
Z1_shift2(:,sz2+1:2*sz2) = Z1_shift1;
Z1_shift2(:,2*sz2+1:3*sz2) = repmat(Z1R,1,sz2);

G1L = G(:,1); G1R = G(:,end);
G2(:,1:sz2) = repmat(G1L,1,sz2);
G2(:,sz2+1:2*sz2) = G;
G2(:,2*sz2+1:3*sz2) = repmat(G1R,1,sz2);

%We do a convolution in the von Mises + log space.
fftZ = fft2(ifftshift(Z1_shift2));
fftG = fft2(ifftshift(G2));
newZZ = fftshift(ifft2(fftZ.*fftG));
newZ = newZZ(:,sz2+1:2*sz2); %cut out original image

%These are the coordinates of the X-Y plane in transformed space.
[T,R] = cart2pol(X,Y);
logR = log(R);

%Now, we interpolate newZ at [T,logR] and remove NaN 
Z2_shift = interp2(newT',newR',newZ',T,logR);
Z2_shift(isnan(Z2_shift))=0;

%Interpolate back to original size input image
[XX,YY] = meshgrid(-cz+1:k_scalar:cz-2);
ZZ = interp2(X,Y,Z2_shift,XX,YY);

%Output is the original size image.
out = ZZ;
