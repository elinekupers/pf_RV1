function [P_svmEnergy, P_svmLinear] = absorptionsClassifierAccuracyStimTemplateWrapper...
                (baseFolder, subFolder, expName, seed)
% Function to classify cone absorption responses with SVM-Energy and SVM-
% Linear template.
% Absorptions are computed by runComputationalObservermodel.m.
%
% Our classifier determines whether the responses belong to a trial with a
% 2-AFC orientation discrimination task (was Gabor clockwise or counter-
% clockwise oriented). This function applies a 2D template to cone responses
% before classifying with a linear SVM.
%
% INPUTS
% baseFolder    :   folder where absorption data live
% subFolder     :   subfolder name when running multiple iterations.
% expName       :   string with experiment name, e.g. 'conedensity' (for
%                   all options see loadExpParams.m)
% seed          :   random number generator seed, in manuscript, we use the
%                   following rule: run1 has rng seed 1, run2 has 2.. etc.%
% Example:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% subFolder  = 'run1';
% expName    = 'conedensity';
% seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
%
% [P_svmEnergy, P_svmLinear] = absorptionsClassifierAccuracyStimTemplateWrapper...
%     (baseFolder, subFolder, expName, seed)

%% 0. Define params
% set random number generator seed and general flags
rng(seed);

% Get experimental params
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
inputType = 'absorptions'; %'absorptions'; % could also be 'current'
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 1:109; % photocurrent responses are temporally delayed
end

% Create folder to save classifier accuracy
saveFolderClassification = fullfile(baseFolder, 'data',  expName, 'classification', inputType, 'stimTemplate', subFolder);
if ~exist('saveFolderClassification', 'dir'); mkdir(saveFolderClassification); end

for eccen = 1:length(expParams.eccentricities)
     
    %% ------------------- Classify absorptions  -------------------
    fprintf('Eccentricity %2.2f\n', expParams.eccentricities(eccen))
    
    % Define filename
    fnameClassify = sprintf(...
        'Classify_coneOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf4.00_lms-%1.1f%1.1f%1.1f.mat',...
        contrasts(end), expParams.polarAngle, sprintf('%i',expParams.eyemovement), expParams.eccentricities(eccen), ...
        expParams.defocusLevels(1), expParams.cparams.noise, expParams.cparams.spatialDensity(2),expParams.cparams.spatialDensity(3),expParams.cparams.spatialDensity(4));
    
    if strcmp(inputType, 'current')
        fnameClassify = ['current_' fnameClassify];
    end
    
    if expParams.verbose
        fprintf('(%s): Classify cone absorption data..\n', mfilename);
        fprintf('(%s): File will be saved as %s\n', mfilename, fnameClassify);
    end
    
    % Preallocate space
    P_svmEnergy  = NaN(1,length(contrasts));
    P_svmLinear  = P_svmEnergy;
    
    for c = 1:length(contrasts)
        fprintf('Contrast %1.4f\n', contrasts(c))
        
        % Load cone responses
        load(fullfile(baseFolder, 'data',  expName, inputType, expName, subFolder, ...
            sprintf('OGconeOutputs_contrast%1.3f_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat', contrasts(c), expParams.eccentricities(eccen))), 'absorptions');
        
        % Get zero contrast stim for mean RGC response
        zeroContrast = load(fullfile(baseFolder, 'data',  expName, inputType, expName, subFolder, ...
            sprintf('OGconeOutputs_contrast0.000_pa0_eye11_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat', expParams.eccentricities(eccen))), 'absorptions');
      
        % Reshape zero contrast stim to a mean noise template
        zeroContrastPermuted = permute(zeroContrast.absorptions,[2 3 4 1 5]);
        zeroContrastReshaped = reshape(zeroContrastPermuted, size(zeroContrastPermuted,1), size(zeroContrastPermuted,2), size(zeroContrastPermuted,3),[]);
        zeroContrastMean = mean(mean(zeroContrastReshaped,3),4);
        
        clear zeroContrast zeroContrastPermuted zeroContrastReshaped
        
        % Truncate time if needed
        if size(absorptions,4) > length(selectTimePoints)
            dataIn = absorptions(:,:,:,selectTimePoints,:);
        else
            dataIn = absorptions;
        end
        
        %% Get template from noiseless stimulus response
        stimTemplate =  getStimTemplateForSVMClassification(baseFolder, subFolder, expName, [], contrasts(c), expParams.eccentricities(eccen), selectTimePoints);
        
        % Classify absorptions
        [P_svmEnergy(c), P_svmLinear(c)] = getClassifierAccuracyStimTemplate(dataIn, stimTemplate, zeroContrastMean);

        if expParams.verbose 
            fprintf('(%s): Classifier accuracy for stim contrast %1.4f is..  Energy: %3.2f\tLinear: %3.2f\n', mfilename, contrasts(c), P_svmEnergy(c), P_svmLinear(c)); 
        end
    end % contrast
    
    % Save
    parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'P_svmEnergy',P_svmEnergy,'P_svmLinear',P_svmLinear, 'expParams', expParams);

end % eccentricities

return