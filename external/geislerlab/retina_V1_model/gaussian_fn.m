function out = gaussian_fn(std_x, std_y, szx, szy)

%gaussian_fn outputs a Gaussian with a peak of 1 and standard deviations of std_x and std_y for a matrix of size
%[szx,szy].

%Place peak at center (defined by "center" function)
[x,y] = meshgrid(-floor(szx/2):ceil(szx/2)-1,-floor(szy/2):ceil(szy/2)-1);

%Prevent division by 0.
if std_x == 0
    std_x = szx*0.001;
end

if std_y == 0
    std_y = szx*0.001;
end

%Output
out = exp(-.5*(x.^2/std_x^2+y.^2/std_y^2));
