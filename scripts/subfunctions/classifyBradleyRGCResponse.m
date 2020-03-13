function P = classifyBradleyRGCResponse(T_pool_stim)

sz=size(T_pool_stim);
dataIn = [];
if length(sz)==4 % [rows, cols, tOrient, nTrials]
    
    % add singleton dimension for time points
    data(1,:,:,:,:) = T_pool_stim; % we have time samples, rows, cols, stim orientations, nTrials
    data = permute(data, [5, 2, 3, 1, 4]); % we want nTrials, rows, cols, time samples, stim orientations
    
elseif length(sz)==3 % [rows, cols, tOrient]
    % add singleton dimension for time points
    
    data(1,1,:,:,:) = T_pool_stim; % we have n trials, time samples, rows, cols, stim orientations
    data = permute(data, [1, 3, 4, 2, 5]); % we want n trials, rows, cols, time samples, stim orientations
else
    error('(%s) input data has too many dimensions!', mfilename)
end

    % Save percent correct
    P = getClassifierAccuracy(data);
    
end


