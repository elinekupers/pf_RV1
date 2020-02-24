function rgcIndices = getRGCMosaic(nrowsCones, ncolsCones, n)
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
%   n   :   linear cone to RGC ratio
%
% OUTPUTS
%   rgc :   rectangular RGC grid with same size as cone mosaic

% % current version needs x,y to be equal (i.e. request square matrix)
% assert(x==y);
% nmax = x^2;

% input cone mosaic
r = 1:n:nrowsCones;
c = 1:n:ncolsCones;

% [c, r] = meshgrid(1:n:ncolsCones, 1:n:nrowsCones);
% rgcIndices = sub2ind([nrowsCones, ncolsCones], r, c);

return

% % debug
% figure; hold all;
% plot(X, Y, 'ko'); hold on;
% plot(X(logical(rgc)), Y(logical(rgc)), 'r.');
% l = findobj(gcf);
% legend([l(end), l(3)], {'Cones', 'RGCs'})
% xlim([0 nrowsCones]); ylim([0 ncolsCones]); axis square; box off;
% title(sprintf('Mosaic for Cones:RGCs = 1:%d', n));
% xlabel('# cones');
% ylabel('# cones');
% set(gca, 'TickDir', 'out', 'FontSize', 12);
