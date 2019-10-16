function out = create_resmap(sx, sy, c_std, ftx, fty)

%create_resmap creates a foveated spatial resolution map of size [sx,sy], measured in pixels. Pixels per degree is fixed
%at 120. It is assumed that the standard deviation of the receptive field center of a foveal ganglion cell is c_std.
%(ftx,fty) represents the coordinates of the fovea, in degrees. The output represents the right eye (nasal to the left).

%Pixels per degree
ppd = 120;

%The resolution map
[X,Y] = meshgrid(1-center(sx):sx-center(sx),sy-center(sy):-1:1-center(sy));
X = X./ppd; Y = Y./ppd;
c_scalar = spacing_fn(0,0)/c_std;
resmap = spacing_fn(X-ftx,Y-fty)./c_scalar; 

%Output
out = resmap;