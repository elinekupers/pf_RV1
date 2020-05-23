function [expName, subFolder, seed] = pfRV1_prepHPC(taskID)

% function to run different experiments with different seeds on NYU's HPC



switch taskID
    
    case {1, 2, 3, 4, 5}
        expName   = 'conedensity';
        subFolder = sprintf('run%d', taskID);
        seed      = taskID;
        
    case {6, 7, 8, 9, 10}
        expName   = 'defaultnophaseshift';
        subFolder = sprintf('run%d', taskID-5);
        seed      = taskID-5;
        
    case {11, 12, 13, 14, 15}
        expName   = 'default';
        subFolder = sprintf('run%d', taskID-10);
        seed      = taskID-10;
        
end

return