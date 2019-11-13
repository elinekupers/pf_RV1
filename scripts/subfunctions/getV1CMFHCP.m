function dataMeridians = getV1CMFHCP()
% Get V1 cortical magnification factor data (CMF, mm2 of cortex / deg
% visual field) between 1 and 6 degrees eccentricity. These data were used
% on Noah Benson's SFN 2019 poster.

%% V1 CMF
T = readtable(fullfile(pfRV1rootPath, 'external', 'data', 'HV-DV-surface-areas.csv'));

allMasks.eccen1_6   = T.min_eccen==1; % deg
allMasks.eccen0_35  = T.min_eccen==0; % deg
allMasks.eccen35_7  = T.min_eccen==3.5; % deg
polang              = T.delta_angle==30; % deg

cmfData = T.surface_area_mm2;

% TO CHECK
temporalIdx = zeros(size(T,1),1);
nasalIdx    = temporalIdx;
inferiorIdx = temporalIdx;
superiorIdx = temporalIdx;

for ii = 1:size(T,1)
    temporalIdx(ii) = (strcmp(T.meridian{ii}, 'horizontal') && strcmp(T.label{ii},'V1d'));
    nasalIdx(ii)    = (strcmp(T.meridian{ii}, 'horizontal') && strcmp(T.label{ii},'V1v'));
    inferiorIdx(ii) = (strcmp(T.meridian{ii}, 'vertical') && strcmp(T.label{ii},'V1d'));
    superiorIdx(ii) = (strcmp(T.meridian{ii}, 'vertical') && strcmp(T.label{ii},'V1v'));
end

% Get integral of V1 CMF for three different eccentricity wedges
fn = fieldnames(allMasks);

dataMeridians = struct();
for jj = 1:3
    
    thisEccenMask = allMasks.(fn{jj});
    
    nasal    = cmfData(nasalIdx & thisEccenMask & polang);
    superior = cmfData(superiorIdx & thisEccenMask & polang);
    temporal = cmfData(temporalIdx & thisEccenMask & polang);
    inferior = cmfData(inferiorIdx & thisEccenMask & polang);
       
    nasalMeanSubject = mean(nasal);
    superiorMeanSubject = mean(superior);
    temporalMeanSubject = mean(temporal);
    inferiorMeanSubject = mean(inferior);

    dataMeridians.(fn{jj}) = [nasalMeanSubject, superiorMeanSubject, temporalMeanSubject, inferiorMeanSubject];

% fprintf('V1 CMF 1-6 deg eccen, 15 deg polar angle wedge')
% fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(dataIn))
% fprintf('Vertical-Meridian Asymmetry (North / South):\t %1.2f%%\n', vma(dataIn))

end


