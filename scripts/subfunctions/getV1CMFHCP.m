function dataMeridians = getV1CMFHCP()
% Get V1 cortical magnification factor data (CMF, mm2 of cortex / deg
% visual field) between 1 and 6 degrees eccentricity. These data were used
% on Noah Benson's SFN 2019 poster.

%% V1 CMF
% T1 = readtable(fullfile(pfRV1rootPath, 'external', 'data', 'distance-ROI-data_angle-stepsize=10.csv'));
T1 = readtable(fullfile(pfRV1rootPath, 'external', 'data', 'distance-ROI-data_angle-stepsize=10_full.csv'));


allMasks.eccen1_6   = (T1.min_eccentricity_deg==1 & T1.max_eccentricity_deg==6); % deg
allMasks.eccen0_35  = T1.min_eccentricity_deg==0; % deg
allMasks.eccen35_7  = T1.min_eccentricity_deg==3.5; % deg
allMasks.eccen05_1  = T1.min_eccentricity_deg==0.5; % deg
allMasks.eccen1_2   = (T1.min_eccentricity_deg==1 & T1.max_eccentricity_deg==2); % deg
allMasks.eccen2_4   = (T1.min_eccentricity_deg==2 & T1.max_eccentricity_deg==4); % deg
allMasks.eccen4_6   = (T1.min_eccentricity_deg==4 & T1.max_eccentricity_deg==6); % deg
allMasks.eccen4_8   = (T1.min_eccentricity_deg==4 & T1.max_eccentricity_deg==8); % deg
allMasks.eccen6_8   = (T1.min_eccentricity_deg==6 & T1.max_eccentricity_deg==8); % deg


polang              = T1.angle_delta_deg==30; % deg

cmfData = T1.surface_area_mm2;

% TO CHECK
horizontalVF = zeros(size(T1,1),1);
verticalVF   = horizontalVF;
upperVF      = horizontalVF;
lowerVF      = horizontalVF;

for ii = 1:size(T1,1)
    horizontalVF(ii) = (strcmp(T1.boundary{ii}, 'horizontal'));
    verticalVF(ii)   = (strcmp(T1.boundary{ii}, 'vertical'));
    upperVF(ii)      = (strcmp(T1.boundary{ii}, 'ventral')); % ventral = upper visual field, corresponding to inferior retina
    lowerVF(ii)      = (strcmp(T1.boundary{ii}, 'dorsal')); % dorsal = lower visual field,  corresponding to superior retina
end

% Get integral of V1 CMF for three different eccentricity wedges
fn = fieldnames(allMasks);

dataMeridians = struct();

for jj = 1:numel(fn)
    
    thisEccenMask = allMasks.(fn{jj});
    
    horz    = cmfData(horizontalVF & thisEccenMask & polang);
    vert = cmfData(verticalVF & thisEccenMask & polang);
    upr = cmfData(upperVF & thisEccenMask & polang);
    lowr = cmfData(lowerVF & thisEccenMask & polang);
       
    dataMeridians.individualSubjects.(fn{jj})   = cat(3, [horz, vert, upr, lowr]);

end


