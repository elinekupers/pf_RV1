% run Performance Field version of Modelfest psychophysics experiment

%% General project path
addpath(genpath(pwd))

% open GL
AssertOpenGL;

% Do you want the debug window?
PsychDebugWindowConfiguration(0, 0.4);
Screen('Preference', 'SkipSyncTests', 1)

% cal  = 'CBI_Propixx'; %
cal  = 'Rm957C_CRT'; %
%cal = 'MacRet10Bit'; % calibration for the new laptop; Note: 60 Hz 

d    = loadDisplayParams(cal);

% Make sure to check which experiment you are running!!
% Go to line 86 in'pfmf_InitStaircaseParams.m' and change the starting stim
% level accordingly!
whichExperiment = 'East'; 

%% Initialize parameters for display, staircase, stimulus, and subject
dataDir          = pfmf_InitDataDir;
display          = pfmf_InitDisplay(cal); 
display          = pfmf_InitFixParams(display);
stimParams       = pfmf_InitStimParams(display, whichExperiment); 

display          = pfmf_InitFixParams(display);
stairParams      = pfmf_InitStaircaseParams(stimParams, display);
subjectParams    = getSubjectParams(dataDir);

trialGenFuncName = 'GaborTrial';

%% Get subject data and log file

logFID(1) = fopen(fullfile(dataDir,[subjectParams.name '.log']), 'at');
fprintf(logFID(1), '%s\n', datestr(now));
fprintf(logFID(1), '%s\n', subjectParams.comment);

if(~isempty(stairParams.curStairVars))
    fprintf(logFID(1), '%s = [ %s ]', stairParams.curStairVars{1}, num2str(stairParams.curStairVars{2}));
end 
fprintf(logFID(1), '\n');
logFID(2) = 1;
hideCursor = false;

%% Do the experiment
devices         = getDevices;
display         = openScreen(display,hideCursor);
display.devices = devices;
priorityLevel   = 0;

results         = doStaircase(display, stairParams, stimParams, trialGenFuncName, priorityLevel, logFID);
display         = closeScreen(display);

save(fullfile(dataDir), 'results')





