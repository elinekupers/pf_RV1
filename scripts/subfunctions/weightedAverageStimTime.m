function dataOut = weightedAverageStimTime(dataIn,interpFilters, meanCur)   
% compute the weighted average over time 
    sz = size(dataIn);
    nTimePoints = size(dataIn,4);
    stim = zeros(1,nTimePoints);
    stim(1:28) = 1;                         % stimulus was on for 28 2-ms time bins (i.e., 56 ms)
    wTime = conv(stim, interpFilters(:,1)); % convolve with l-cone temporal filter
    wTime = wTime(1:nTimePoints);           % limit to length of simulated data
    wTime = wTime / sum(wTime);             % normalize to sum of 1 for weighted average
    
    % initialize an array to store the weighted sum of the current data for
    % this run
    dataOut = zeros(sz(1), sz(2), sz(3), 1, sz(5));
    
    % computed the weighted average
    for ii = 1:nTimePoints
        dataOut = dataOut + dataIn(:,:,:,ii,:)*wTime(ii);
    end
    
    dataOut = dataOut ./sum(dataOut(:));
    dataOut = dataOut .* meanCur(1);
    
    return