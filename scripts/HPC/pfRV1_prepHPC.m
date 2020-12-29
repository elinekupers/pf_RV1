function [expName, subFolder, seed] = pfRV1_prepHPC(taskID,varargin)

% function to run different experiments with different seeds on NYUs HPC
if nargin<2
    prefixSubfolder = '';
else
    prefixSubfolder = varargin{1};
end



switch taskID
    
    case {1, 2, 3, 4, 5}
        expName   = 'conedensity';
        subFolder = sprintf('run%d%s', taskID,prefixSubfolder);
        seed      = taskID;

    case {6, 7, 8, 9, 10}
        expName   = 'defaultnophaseshift';
        subFolder = sprintf('run%d%s', taskID-5, prefixSubfolder);
        seed      = taskID-5;
        
    case {11, 12, 13, 14, 15}
        expName   = 'default';
        subFolder = sprintf('run%d%s', taskID-10, prefixSubfolder);
        seed      = taskID-10;
        
end

return