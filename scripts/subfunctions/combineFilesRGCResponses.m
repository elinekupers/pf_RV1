function combineFilesRGCResponses(baseDir, run)
% Function to combine matfiles with simulated RGC responses that were too 
% big (7.49 GB) to upload to OSF storage as is (max 5.2 GB). 
% To make sure all data are accessible, we separated the first dimension of
% the rgcResponse arrays (nr trials) into half: trials 1-50 and 50-100.
% This function will combine the two halves back into one big file.
% 
% Note: The files are only too big for ratio a 2:1 cone:RGC ratio, at the 
% fovea (eccen = 0.00 deg).
%
%
% Example:
%   baseDir = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/data/conedensity/rgc/';
%   run     = 1;
%   combineFilesRGCResponses(baseDir, run)


%% Specify exp params
ratio     = 1;
eccen     = 0.00;
contrasts = [0, 0.0050, 0.0100, 0.0150, 0.0200, 0.0250, 0.0300, 0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000];

for c = 1:length(contrasts)
    
    % define name
    fname = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%1.2f_absorptionrate', ratio, contrasts(c), eccen);
    
    % load name
    load(fullfile(baseDir, sprintf('run%d/ratio%d/',run, ratio), [fname '_trial1-50.mat']));
    load(fullfile(baseDir, sprintf('run%d/ratio%d/',run, ratio), [fname '_trial51-100.mat']));
    
    % separate first 50 and last 50 trials (first dim)
    rgcResponse  = cat(1, rgcResponse1_50, rgcResponse51_100);
    
    save(fullfile(baseDir, sprintf('run%d/ratio%d/',run, ratio), [fname '.mat']), 'rgcResponse', 'expParams', 'rgcParams', 'contrasts', '-v7.3')

end

return