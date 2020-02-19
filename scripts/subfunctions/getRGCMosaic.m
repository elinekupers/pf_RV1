function rgc = getRGCMosaic(x, y, n)
%
% Function to build rectangular RGC mosaic by subsampling from cone mosaic.
% Assumptions:
%   - center RGC overlaps with center cone
%   - cones and RGC have rectangular spacing
%   - cone and RGC mosaic are square matrices
%
% INPUTS
%   x   :   integer of cones on x-axis
%   y   :   integer of cones on y-axis
%   n   :   cone to RGC ratio
%
% OUTPUTS
%   rgc :   rectangular RGC grid with same size as cone mosaic

% current version needs x,y to be equal (i.e. request square matrix)
assert(x==y);
nmax = x^2;

% input cone mosaic
[X,Y] = meshgrid(1:x, 1:y);

if n == 1
    % create empty rgc matrix with the same sample rate as cone matrix
    rgc = ones(x,y);
    
elseif n == 2
    
    % make checkerboard for ratio of 2
    rgc = zeros(x,y);
    
    % Change phase offset (start of indices) depending if even/uneven
    % matrix
    if mod(x,2)~=0
        indices = 2:n:nmax;
    else
        indices = [2:n:x; (x+1):n:(2*x)];
        indicesAll = [];
        
        for ii = 0:((y/2)-1)
            indicesAll  = [indicesAll; indices+(2*x*ii)];
        end
        rgc(indicesAll) = 1;
    end
    
    rgc(indices)=1;
    
else
     % downsample grid in some long-windy way
    rgc = zeros(x,y);
    
    if mod(log2(n),1) == 0
        n_sub = log2(n);
    else
        n_sub = ceil(n/2);
    end
    
    if mod(x,2)==0
        rx = X(1,1:n_sub:end);
        ry = Y(1:n_sub:end,1);
    elseif mod(x,2)~=0
        n_sub = log2(n);
        rx = X(1,1:n_sub:end);
        ry = Y(2:n_sub:end,1);
    end
    rgc(rx,ry) = 1;
end


% debug
figure; hold all;
plot(X, Y, 'ko'); hold on;
plot(X(logical(rgc)), Y(logical(rgc)), 'r.');
l = findobj(gcf);
legend([l(end), l(3)], {'Cones', 'RGCs'})
xlim([0 x]); ylim([0 y]); axis square; box off;
title(sprintf('Mosaic for Cones:RGCs = 1:%d', n));
xlabel('# cones');
ylabel('# cones');
set(gca, 'TickDir', 'out', 'FontSize', 12);

return
