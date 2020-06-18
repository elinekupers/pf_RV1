function CMF = HortonHoytCMF(eccentricity)
% CMF as a function of eccentricity, from Horton & Hoyt (1991). The
% representation of the visual field in human striate cortex. A revision of
% the classic Holmes map. Arch Ophthalmol. 1991 Jun;109(6):816-24.

% CMF in mm2/deg2, eccentricity in deg

%     CMF_mmperdeg = 17.3./(eccentricity+0.75);
%     CMF = CMF_mmperdeg.^2;
    
     CMF = 300./((eccentricity+0.75).^2);

end