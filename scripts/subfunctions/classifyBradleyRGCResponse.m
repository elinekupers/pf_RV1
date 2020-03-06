function P = classifyBradleyRGCResponse(T_pool_stim, saveDir, saveData)

[rows, cols, tLoc, tOrient, tEccen, nTrials] =size(T_pool_stim);

% we need dims: trials x rows x cols x time samples x stimuli) 
data = permute(T_pool_stim, [6, 1, 2, 3, 5, 4]);

for te = 1:tEccen
    
    tmpData = squeeze(data(:,:,:,:,te,:));
    
    for tl = 1:tLoc
        
        dataIn = tmpData(:,:,:,tl,:);
        
        % Save percent correct
        P(te,tl) = getClassifierAccuracy(dataIn);
        
    end
    
    if saveData
        if ~exist(fullfile(saveDir, 'classification', 'rgcBradley'), 'dir')
            mkdir(fullfile(saveDir, 'classification', 'rgcBradley'));
        end
        save(fullfile(saveDir, 'classification', 'rgcBradley', 'classify_rgcResponse_Bradley.mat'), 'P', '-v7.3')
    end
    
end