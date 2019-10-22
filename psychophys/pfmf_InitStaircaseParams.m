function stairParams = pfmf_InitStaircaseParams(stimParams, display)

% name of expt
stairParams.conditionName        =  {'contrast detect'}; 

% how many independent staircases
numStairs = 1;

% these correspond to the different staircases 
stairParams.stimDuration          = stimParams.duration;

% % how often to pause?
% stairParams.pause_between_blocks  = 30; % pause every 30 trials

% decision subject will make
stairParams.alternativeVarName   = 'present_or_absent';

% decision values
stairParams.alternativeVarValues  = [1 0]; % 1 is present, 0 is absent

% decision keys
stairParams.responseSet          = 'as';

% variable that is adjusted by staircase
stairParams.adjustableVarName   =  'contrast';    
stairParams.adjustableVarValues =  logspace(0,-3,30); %[0 logspace(-4,-1.4,29)];

% put things in here to run the staircase separately for a given condition
stairParams.curStairVars        = {}; %{'testPosition', 0}; %[0 90 180 270]}; % polar angle positions (0, 90, 180, 270) in deg

% VARIABLE PARAMETERS
stairParams.randomVars = {}; 

% limit expt in case of lack of convergence
stairParams.maxNumTrials            = 120;

% end expt after this many reversals
stairParams.maxNumReversals         = 120;

% increment size for correct answers for each successive reversal (normally
% these numbers should go down as num reversals incr)
stairParams.numCorrectForStep       = 1;
stairParams.numIncorrectForStep     = 1;

stairParams.incorrectStepSize       = -1;
stairParams.correctStepSize         = -1;

% auditory feedback?
stairParams.feedback                = 'auditory'; %{'click' 'none')
stairParams.responseSet             = 'as'; 

% display dur of stimulus on control screen
stairParams.showTiming              = false;

% intertrial interval in seconds
stairParams.iti  = 2;



%%
% This specifies the intitial value of the staircase, as an index into
% stairParams.alternativeVarValues. If there are multiple interleaved
% staircases, then we separately set the intial value for each staircase.

initIndex = repmat(length(stairParams.adjustableVarValues), size(stairParams.curStairVars{2})); 

stairParams.adjustableVarStart = initIndex; 


end