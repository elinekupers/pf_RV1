function colors = pfRV1_getColors(nrColors)

% This function will return a matrix with either 4 RGB color codes 
% (red, blue, green, black) or 5 RGB colors codes (red, blue, green, gray, black).

switch nrColors
    case 4
        colors = [166,206,227; 31,120,180; 178,223,138; 51,160,44]/255;
    case 5
        colors = [0, 0, 4; 129, 37, 129; 229, 89, 100; 181, 54, 122; 251, 135, 97]/255;
        
end