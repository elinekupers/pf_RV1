function CMF = RovamuVirsuCMF(eccentricity)
% CMF for the cardinal merdians, from Romavo & Virsu (1979).
% An estimation and application of the human cortical magnification factor.
% Exp. Brain. Rs. vol. 37, pp. 495-510.
%
% CMF in mm2/deg2 for visual field meridians, eccentricity in deg.

% CMF at fovea
M0 = 7.99; % mm/deg

% CMF functions for the different cardinal meridians (in mm/deg)
nasal = M0./(1+ (0.33 .* eccentricity) + (0.00007.*eccentricity.^3));
superior = M0./(1+ (0.42 .* eccentricity) + (0.00012.*eccentricity.^3));
temporal = M0./(1+ (0.29 .* eccentricity) + (0.000012.*eccentricity.^3));
inferior = M0./(1+ (0.42 .* eccentricity) + (0.000055.*eccentricity.^3));

CMF.nasalVF = nasal.^2; % mm2/deg2
CMF.superiorVF = superior.^2; % mm2/deg2
CMF.temporalVF = temporal.^2; % mm2/deg2
CMF.inferiorVF = inferior.^2; % mm2/deg2

CMF.nasalR = temporal.^2; % mm2/deg2
CMF.superiorR = inferior.^2; % mm2/deg2
CMF.temporalR = nasal.^2; % mm2/deg2
CMF.inferiorR = superior.^2; % mm2/deg2
CMF.eccentricity = eccentricity; % deg
CMF.comment = 'eccentricity in deg. CMF in mm2/deg2';

end