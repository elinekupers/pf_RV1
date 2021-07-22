function data4D = pf5Dto4D(data5D)
%
% Inputs are 5D: 
%   trials x cols x rows x time points x stimulus types
% Outputs are 4D, folding stimulus types into trials:
%    cols x rows x time points x trials
%
% We might also need pf4Dto5D??

cRows   = size(data5D, 2);
cCols   = size(data5D, 3);
cTime   = size(data5D, 4);
permuteddata5D = permute(data5D, [2, 3, 4, 1, 5]);
data4D = reshape(permuteddata5D, cRows, cCols, cTime, []);

end