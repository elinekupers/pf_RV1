function dataMeridians = getV1CMFHCP(wedgeWidth)
% Function to compute V1 cortical magnification factor data (CMF) in mm2
% cortex/deg2 visual field) between 1 and 6 degrees eccentricity.
%
% These data come from the published article:
%   Benson, Jamison, Arcaro, Vu,Glasser, Coalson, VanEssen, Yacoub, Ugurbil,
%   Winawer, Kay. 2018 The Human Connectome Project 7 Tesla retinotopy
%   dataset: Description and population receptive field analysis.
%   Journal of Vision, 18:23. doi:10.1167/18.12.13.
% 
% Written by EK 2019 NYU

%% 0. Read data in
T1 = readtable(fullfile(pfRV1rootPath, 'external','data', 'DROI_table.csv'));
cmfData = T1.surface_area_mm2;

numSubjects = size(T1,1);

%% 1. Define eccentricity, polar angle data masks
allMasks.eccen1_2   = (T1.min_eccentricity_deg==1 & T1.max_eccentricity_deg==2); % deg
allMasks.eccen2_3   = (T1.min_eccentricity_deg==2 & T1.max_eccentricity_deg==3); % deg
allMasks.eccen3_4   = (T1.min_eccentricity_deg==3 & T1.max_eccentricity_deg==4); % deg
allMasks.eccen4_5   = (T1.min_eccentricity_deg==4 & T1.max_eccentricity_deg==5); % deg
allMasks.eccen5_6   = (T1.min_eccentricity_deg==5 & T1.max_eccentricity_deg==6); % deg

polang              = T1.angle_delta_deg==wedgeWidth; % deg

%% 2. Get CMF data per polar angle

% preallocate space
horizontalVF = NaN(numSubjects,1);
verticalVF   = horizontalVF;
upperVF      = horizontalVF;
lowerVF      = horizontalVF;

% Get indices for visual fields
for ii = 1:numSubjects
    horizontalVF(ii) = (strcmp(T1.boundary{ii}, 'horizontal'));
    verticalVF(ii)   = (strcmp(T1.boundary{ii}, 'vertical'));
    upperVF(ii)      = (strcmp(T1.boundary{ii}, 'ventral')); % ventral = upper visual field, corresponding to inferior retina
    lowerVF(ii)      = (strcmp(T1.boundary{ii}, 'dorsal')); % dorsal = lower visual field,  corresponding to superior retina
end

% Use indices to get data for visual fields
fn = fieldnames(allMasks);

% Get data struct
dataMeridians = struct();

for jj = 1:numel(fn)
    % Select eccentricity mask
    thisEccenMask = allMasks.(fn{jj});
    % Select corresponding CMF data for visual field, eccen, polar ang
    horz    = cmfData(horizontalVF & thisEccenMask & polang);
    vert    = cmfData(verticalVF & thisEccenMask & polang);
    upr     = cmfData(upperVF & thisEccenMask & polang);
    lowr    = cmfData(lowerVF & thisEccenMask & polang);
    
    % Concatenate data
    dataMeridians.individualSubjects.(fn{jj})   = cat(3, [horz, vert, upr, lowr]);
end

return
