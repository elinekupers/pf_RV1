function asymmetryPercent = hva(dataIn)
% Function to calculate horizontal-vertical asymmetry. 
%
% INPUT:
%   dataIn            : (vector) meridian data. Data have contain cardinal
%                           meridians on the retina, in the order:
%                           (1) EAST (nasal/temporal retina, depending on eye)
%                           (2) NORTH (superior retina)
%                           (3) WEST (nasal/temporal retina, depending on eye)
%                           (4) SOUTH (inferior retina)
%                           
%                        or cardinal and off-cardinals in the order:
%                           (1) EAST  (nasal/temporal retina, depending on eye)
%                           (2) NORTHEAST 
%                           (3) NORTH (superior retina) 
%                           (4) NORTHWEST,
%                           (5) WEST  (nasal/temporal retina, depending on eye) 
%                           (6) SOUTHWEST, 
%                           (7) SOUTH (inferior retina)
%                           (8) SOUTHEAST
%
% OUTPUT:
%   asymmetryPercent  : horizontal-vertical asymmetry in percent difference

if numel(dataIn) == 4  % If data contain cardinals only
    horz = [1 3];
    vert = [2 4];
    
elseif numel(dataIn) == 8 % If data contain cardinals and off cardinals
    horz = [1 5];
    vert = [3 7];
    
else
    error('%s, number of elements in data (%d) do not have correct size', mfilename, numel(sz))
end

asymmetry = (mean(dataIn(horz)) - mean(dataIn(vert))) ./ mean(dataIn([horz, vert]));
asymmetryPercent = 100*asymmetry;


return