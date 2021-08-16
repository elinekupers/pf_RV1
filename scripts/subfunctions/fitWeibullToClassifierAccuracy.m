function weibullFit = fitWeibullToClassifierAccuracy(contrastLevels, accuracy, ...
                        initWeibullParams, xUnits)
% INPUTS:
% contrastLevels      :     row vector with Michelson stimulus contrast 
%                           (units are a fraction of 1) where 0 is 0% (no 
%                           contrast) and 1 is 100% contrast.
% accuracy            :     column vector with classifier accuracy per  
%                           stimulus contrast level, in % correct, ranging
%                           from chance (~50% for a 2-AFC task) to 100%. 
% initWeibullParams   :     inital slope and threshold for first stage fitting
%                           example would be [3 0.01]
%                           beta = 3, predefined in ogFitWeibull, gives 
%                           threshold of Weibull fit at ~80% correct.
% xUnits              :     row vector with finely sampled x-axis to fit 
%                           final weibull, this allows for a more accurate
%                           estimation of the contrast threshold.
%
% OUTPUTS
% weibullFit          :     struct with the following params related to the
%                           weibull that has been fitted and where "ctr"
%                           stands for contrast:
%                           - init:    stored initial slope and threshold
%                           - ctrvar:  fitted contrast variables, i.e.,
%                             slope and threshold from first stage fit.
%                           - ctrpred: full weibull prediction with finely
%                             sampled x-axis. 
%                           - ctrthresh: contrast threshold at ~80% (0.5^1/3) 
%                             in stimulus contrast using ctrpred.

weibullFit = struct();
weibullFit.init = initWeibullParams;
weibullFit.ctrvar    = [];
weibullFit.ctrpred   = [];
weibullFit.ctrthresh = [];
weibullFit.data      = [];

% Make a Weibull function first with contrast levels and then search for
% the best fit with the classifier data
weibullFit.ctrvar = fminsearch(@(x) ogFitWeibull(x, contrastLevels, ...
                        accuracy, 100), weibullFit.init);

% Then fit a Weibull function again, but now with the best fit parameters
% from the previous step, and finely sampled x-axis
weibullFit.ctrpred = 100*ogWeibull(weibullFit.ctrvar, xUnits); % multiply by 100 to get percent correct

% Extract contrast threshold, and save corresponding data
weibullFit.ctrthresh = weibullFit.ctrvar(2);
weibullFit.data  = accuracy;

end