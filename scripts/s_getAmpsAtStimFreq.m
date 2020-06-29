%% s_getAmpsAtStimFreq.m

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
subFolder  = 'run5';

expName   = 'conedensity';
expParams = loadExpParams(expName, false);
contrasts = expParams.contrastLevels;

% Define DoG Params
rgcParams.DoG.kc     = 1/3;                 % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;                % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc;  % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = 'absorptionrate';   % are we dealing with cone absorptions or current?
rgcParams.fov        = 2;                   % dva
rgcParams.stimSF     = expParams.spatFreq;  % cpd
selectTimePoints     = 1:28;                % timepoints (2 ms samples)
sfOfStim             = rgcParams.stimSF * rgcParams.fov;

cone2RGCRatios       = 1:5;
nrEccen = length(expParams.eccentricities);

for ii = 1:length(cone2RGCRatios)
    
    % Subsampling ratio
    rgcParams.cone2RGCRatio = cone2RGCRatios(ii);
    
    % Center Gauss RGC
    sigma.center = rgcParams.cone2RGCRatio*rgcParams.DoG.kc;
    
    % Surround Gauss RGC
    sigma.surround = rgcParams.cone2RGCRatio*rgcParams.DoG.ks;
    
    for eccen = 1:nrEccen
        
        for c = 1:length(contrasts)
            load(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ii), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%1.2f_%s.mat', ii,  contrasts(c), expParams.eccentricities(eccen), rgcParams.inputType)), 'rgcResponse');
            
            data = rgcResponse(:,:,:,selectTimePoints,:);
            
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
            
            % take sum across time 
            % makes a 3D matrix: rows x cols x (trials x stimulus)
            data = sum(data,4);
            
            % take difference between cw and ccw
            dataDiff = data(:,:,1:nTrials) - data(:,:,(nTrials+1):end);
            dataDiffMean = mean(dataDiff,3);
            
            % Compute fourier transform the cone array outputs
            amps  = abs(fft2(dataDiffMean));
            
            % Get amps at stim freq
            x      = (0:nrows-1)/nrows*rgcParams.fov; % degrees
            fs     = (0:nrows-1)/2;
            stimIdx = (nrows+1)-sfOfStim;
            ampsStim(c,eccen,ii) = amps(3,stimIdx);
        end
    end
    
end


for ii = 1:length(cone2RGCRatios)
     for eccen = 1:nrEccen
         [a(ii,eccen), ai(ii,eccen)] = max(squeeze(ampsStim(:,eccen, ii)));
     end
end

for ii = 1:length(cone2RGCRatios) 
    [b(ii), bi(ii)] = max(a(ii,:));
end

dataFolder = fullfile(baseFolder, 'data',  expName, 'amplitudes');
if ~exist('dataFolder', 'dir'); mkdir(dataFolder); end
save(fullfile(dataFolder, sprintf('amplitudesAtStim_ratio1_5_%s', subFolder)), 'data', 'amps', 'ampsStim');


