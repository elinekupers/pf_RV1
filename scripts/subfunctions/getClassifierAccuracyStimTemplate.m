function [PEnergy, PLinear] = getClassifierAccuracyStimTemplate(data, stimTemplate, zeroContrastMean)
% Function to train and test linear SVM classifier on cone data with 
% cross-validation. We apply a noiseless stimulus template to the responses
% before applying the 2D FFT, to see if it improves classifier performance.
%
%       [PEnergy, PLinear] = getClassifierAccuracyStimTemplate(data)
%
% INPUTS: 
%   data             : 5 dimensional array (trials x rows x cols x time x
%                       x stimuli) 
%   stimTemplate     : struct with noiseless cone/rgc response (rows x cols)
%                       in quadriture phases
%   zeroContrastMean : single frame of only noise (zero contrast) RGC
%                       responses averaged across trials and time (rows x cols)
% OUTPUTS:
%   PEnergy          : classifier accuracy of template-energy observer model in
%                   percent correct for given absorption dataset
%   PLinear          : classifier accuracy of template-linear observer model in
%                   percent correct for given absorption dataset


% Get dimensions of data
fprintf('(%s): Loading and classifying\n', mfilename);
% Get the trials and samples (should be the data for all data sets though
nStimuli = size(data,5);
nTrials  = size(data,1) * nStimuli/2;
tSamples = size(data,4);
nrows    = size(data,2);
ncols    = size(data,3);

%   permute to trials x stimuli x rows x cols x time points
data = permute(data, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
data = reshape(data, [], nrows, ncols, tSamples);

% permute to rows x cols x (trials x stimuli) x time points
data  = permute(data, [2 3 1 4]);

% Remove average
data = data - zeroContrastMean;

% Reshape to (cols x rows) x (trials x stimul) x time points
data = reshape(data, [], size(data,3), size(data,4));

% Apply Energy or Linear template 
for nn = size(data,2):-1:1 % trials
    for tt = size(data,3):-1:1 % time
        x = data(:,nn,tt)';
        
        energyXform(nn,tt) = ...
            (x * stimTemplate.cw_ph1(:))^2 + ...
            (x * stimTemplate.cw_ph2(:))^2 - ...
            (x * stimTemplate.ccw_ph1(:))^2 - ...
            (x * stimTemplate.ccw_ph2(:))^2; 
        
        linearXform(nn,tt) = ...
            x * (stimTemplate.cw_ph2(:) - stimTemplate.ccw_ph2(:));
    end
end

% permute the trial order within each of the two classes
idx = [randperm(nTrials) randperm(nTrials)+nTrials];
energyXform = energyXform(idx, :);
linearXform = linearXform(idx, :);

% Create labels for SVM linear classifier
label = [ones(nTrials, 1); -ones(nTrials, 1)];

% Fit the SVM model.
cvmdlEnergy = fitcsvm(energyXform, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
cvmdlLinear = fitcsvm(linearXform, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);

% predict the data not in the training set.
classLossEnergy = kfoldLoss(cvmdlEnergy);
classLossLinear = kfoldLoss(cvmdlLinear);

% Get percent accuracy
PEnergy = (1-classLossEnergy) * 100;
PLinear = (1-classLossLinear) * 100;

% visualize beta's
% betas = reshape(cvmdlEnergy.Trained{1}.Beta, [nrows, ncols, tSamples]);
% mn_betas = squeeze(mean(betas,3));
% imagesc(fftshift(mn_betas)); box off; set(gca, 'TickDir', 'out', 'FontSize',12); 
% colormap gray; axis image; colorbar;
% 
% title(sprintf('FFT at input freq: %1.3f x10^6', mn_betas(8,3)*10^6));
% set(gca,'CLim', 1.0e-03 .*[-0.3624,0.3624]);
%     savePth = fullfile(pfRV1rootPath, 'figures');
%     print(fullfile(savePth, 'classifierWeightsSVMEnergy_averagedAcrossTime.eps'),'-depsc')

       
return