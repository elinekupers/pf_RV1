function [] = RGCmodel_filterAndClassify()
% This function loads cone absorptions and current for the no-eye
% movement L-only cone mosaic. It then applies our late noise RGC model:
% First filter with DOG, then add white noise, then downsample.
% Finally, we run a linear SVM classifier to get 2-AFC (CW/CCW) accuracy 
% for each data type (absorptions, current, filtered, late noise,
% downsampled 1-5)
%
% See s_2dfilterAndClassify and s_1dfilterAndClassify
runnum = 2;

% pth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/';
pth = fullfile(ogRootPath, 'data');
expname = 'conedensitynophaseshiftlonly500'; %'defaultnophaseshiftlonly500';
subfolder = sprintf('run%d', runnum);
expParams = loadExpParams(expname);

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = false; % print figures or not
rgcParams.saveFigs   = false;
rgcParams.expName    = expname;
rgcParams.subFolder  = [];

% Define DoG Params
ratio = 1;
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;               % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = 'current';          % are we dealing with cone absorptions or current?
rgcParams.cone2RGCRatio = ratio;           % linear ratio
rgcParams.seed       = runnum;

% Define late noise level and downsample factors
lateNoiseLevel   = 1; % std
downSampleFactor = [1 2 3 4 5];    % downsample factor
datatypes = {'absorptions', 'current', 'Filtered', 'LateNoise','DownSampled1',...
    'DownSampled2','DownSampled3' 'DownSampled4' 'DownSampled5'};

hpcArrayID = str2double(getenv('SLUM_ARRAY_TASK_ID'));
if ~isnan(hpcArrayID)
    expParams = hpcArrayID2eccen(hpcArrayID, expParams);
end
for ec = 1:length(expParams.eccentricities)
    for c = 1:length(expParams.contrastLevels)
        contrast = expParams.contrastLevels(c);
        fprintf('Contrast %1.4f..\n',contrast)
        
        expstr = sprintf('contrast%1.4f_pa0_eye00_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat', ...
            contrast, expParams.eccentricities(ec));
        
        % load the simulated photocurrent data
        load(fullfile(pth, 'conecurrent', expname, subfolder, ...
            sprintf('currentMn_OGconeOutputs_%s', expstr)));
        
        % load the simulated absorption data
        load(fullfile(pth, 'coneabsorptions', expname, subfolder, ...
            sprintf('Mn_OGconeOutputs_%s', expstr)));
            
            %% filter the data
            x.absorptions = pf5Dto4D(absorptions);
            x.current     = pf5Dto4D(current);
            clear absorptions current;
            
            [x.rgcResponse, ~, ~, x.Filtered] = rgcLayerLinear(x.current, rgcParams, expParams);
            
            %% next we need to
            % - classify the filtered data
            % - then add noise
            % - then subsample at a few rates and classify again for each subsampling
            % then if it makes sense, run with multiple contrasts on multiple runs
            [rows, cols, nTimePoints, numTrials] = size(x.Filtered);
            
            latenoise        = randn(rows, cols,nTimePoints, numTrials)* lateNoiseLevel;
            x.LateNoise      = x.Filtered + latenoise;
            
            for ii = 1:length(downSampleFactor)
                
                coneArray = zeros(rows, cols);
                rowIndices = 1:downSampleFactor(ii):rows;
                colIndices = 1:downSampleFactor(ii):cols;
                
                x.(sprintf('DownSampled%d',ii)) = x.LateNoise(rowIndices, colIndices,:,:); % cone to RGC sampling + late noise
            end
            
            %     figure(99); clf; hold all;
            for jj = 1:length(datatypes)
                fprintf('%d..',jj)
                data = x.(datatypes{jj});
                data = permute(data, [4, 1, 2, 3]);
                data = reshape(data, numTrials, []);
                labels = [ones(numTrials/2,1);zeros(numTrials/2,1)];
                cvmdl = fitcsvm(data, labels, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
                classLoss = kfoldLoss(cvmdl);
                PercentCorrect(c,jj) = (1-classLoss) * 100;
                
                %     bar(jj,PercentCorrect(jj)'); set(gca, 'FontSize', 12)
                %     ylim([50 100])
                
            end
            fprintf('Done!\n')
            
            % save classifier accuracy data
            saveStr = sprintf('rgcResponses_latenoiselevel%1.1f_withDownsampling_c%d_eccen%1.2f.mat',lateNoiseLevel,c,expParams.eccentricities(ec));
            accuracy = PercentCorrect(c,:);
            save(fullfile(pth,'conecurrent', expname, subfolder, saveStr), 'x','accuracy', 'expParams', 'rgcParams')
            clear accuracy
    end
    
    % save classifier accuracy data
    saveStr = sprintf('classifierAccuracy_latenoiselevel%1.1f_eccen%1.2f_withDownsampling.mat', ...
        lateNoiseLevel, expParams.eccentricities(ec));
    save(fullfile(pth,'conecurrent', expname, subfolder, saveStr), 'PercentCorrect', 'expParams', 'rgcParams')
   
end

return
