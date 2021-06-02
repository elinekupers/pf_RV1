function dataMeridians = getV1SurfaceAreaHCP(wedgeWidth)
% Function to compute V1 surface area in mm2 on the cortex. 
% Note: Not yet CMF! Divide by deg2 visual field for each eccentricity wedge
% between 1 and 6 degrees eccentricity.
%
% These data come from the published article:
%   Benson, Jamison, Arcaro, Vu,Glasser, Coalson, VanEssen, Yacoub, Ugurbil,
%   Winawer, Kay. 2018 The Human Connectome Project 7 Tesla retinotopy
%   dataset: Description and population receptive field analysis.
%   Journal of Vision, 18:23. doi:10.1167/18.12.13.
% 
% Written by EK 2019 NYU

%% 0. Read data in
T1 = readtable(fullfile(pfRV1rootPath, 'external','data', 'benson2021','DROI_table.csv'));
surfaceArea = T1.surface_area_mm2;
subjIDs     = unique(T1.sid);
numSubjects = length(subjIDs);
%% 1. Define eccentricity, polar angle data masks
allMasks.eccen1_2   = (T1.min_eccentricity_deg==1 & T1.max_eccentricity_deg==2); % deg
allMasks.eccen2_3   = (T1.min_eccentricity_deg==2 & T1.max_eccentricity_deg==3); % deg
allMasks.eccen3_4   = (T1.min_eccentricity_deg==3 & T1.max_eccentricity_deg==4); % deg
allMasks.eccen4_5   = (T1.min_eccentricity_deg==4 & T1.max_eccentricity_deg==5); % deg
allMasks.eccen5_6   = (T1.min_eccentricity_deg==5 & T1.max_eccentricity_deg==6); % deg

polang              = T1.angle_delta_deg==wedgeWidth; % deg
hemis(:,1)          = strcmp(T1.hemisphere, 'lh');
hemis(:,2)          = strcmp(T1.hemisphere, 'rh');

%% 2. Get CMF data per polar angle

% preallocate space
horizontalVF = NaN(size(T1,1),numSubjects);
verticalVF   = horizontalVF;
upperVF      = horizontalVF;
lowerVF      = horizontalVF;
subject      = horizontalVF;

% Get indices for visual fields
for ii = 1:numSubjects
    subject(:,ii)      = (T1.sid == subjIDs(ii));
    
    horizontalVF(:,ii) = subject(:,ii) & strcmp(T1.boundary, 'horizontal');
    verticalVF(:,ii)   = subject(:,ii) & strcmp(T1.boundary, 'vertical');
    upperVF(:,ii)      = subject(:,ii) & strcmp(T1.boundary, 'ventral'); % ventral = upper visual field, corresponding to inferior retina
    lowerVF(:,ii)      = subject(:,ii) & strcmp(T1.boundary, 'dorsal'); % dorsal = lower visual field,  corresponding to superior retina
end

% Use indices to get data for visual fields
fn = fieldnames(allMasks);

% Get data struct
dataMeridians = struct();

for jj = 1:numel(fn)
    % Select eccentricity mask
    thisEccenMask = allMasks.(fn{jj});
    
    for kk = 1:numSubjects
        for h = 1:size(hemis,2)
            % Select corresponding CMF data for visual field, eccen, polar ang
            horz(:,kk,h)    = surfaceArea(horizontalVF(:,kk) & thisEccenMask & polang & hemis(:,h));
            vert(:,kk,h)    = surfaceArea(verticalVF(:,kk) & thisEccenMask & polang & hemis(:,h));
            upr(:,kk,h)     = surfaceArea(upperVF(:,kk) & thisEccenMask & polang & hemis(:,h));
            lowr(:,kk,h)    = surfaceArea(lowerVF(:,kk) & thisEccenMask & polang & hemis(:,h));
        end
    end

    % Concatenate data (4 cardinal meridians, 163 subjects, 2 hemi (lh/rh)
    dataMeridians.individualSubjects.(fn{jj})   = cat(1, horz, vert, upr, lowr);
end

return
