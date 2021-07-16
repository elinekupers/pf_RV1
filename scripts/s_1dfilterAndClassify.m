numTimePoints   = 64;  
numTrials       = 400;  % Number of trials per condition
numExperiments  = 2;    % number of experiments to loop over
earlyNoiseLevel = 2;    % (std)
lateNoiseLevel  = 1;    % 
f               = 2;    % Signal frequency
downSampleFactor= [2 4 6];    % downsample factor
nonlins         = {'ThresholdAndSquare' 'Normalization' 'none'};
which_nonlin    = nonlins{3};

t       = (1:numTimePoints)';
signal1 = sin(t*2*pi/numTimePoints*f); 
signal2 = cos(t*2*pi/numTimePoints*f);
s       = [signal1 * ones(1,numTrials) signal2 * ones(1,numTrials)];
labels  = [ones(numTrials,1); zeros(numTrials,1)];

figure(1), clf; 

 
PercentCorrect = NaN(5,numExperiments);
x = [];

for ii = 1:numExperiments
    earlynoise     = randn(numTimePoints,numTrials*2)*earlyNoiseLevel;
    latenoise      = randn(numTimePoints,numTrials*2)* lateNoiseLevel;

    x.Template     = (signal1-signal2)' * (s+earlynoise);
    x.Noisy        = s+earlynoise; % analagous to noisy photon absorptions 
    x.Filtered     = bandpass(x.Noisy,f .* [.5 2],numTimePoints); % analagous to DoG spatial filter of RGCs 

    x.NonLinearity = applynonlinearity(x.Filtered, which_nonlin); % analagous to RGC nonlinearity (eg membrane current to spiking)

    %x.LateNoise    = x.NonLinearity + latenoise;

    x.DownSampled1 = downsample(x.NonLinearity+latenoise, downSampleFactor(1)); % cone to RGC sampling + late noise
    x.DownSampled2 = downsample(x.NonLinearity+latenoise, downSampleFactor(2)); % cone to RGC sampling + late noise
    x.DownSampled3 = downsample(x.NonLinearity+latenoise, downSampleFactor(3)); % cone to RGC sampling + late noise
    
    tt.Template     = 0;
    tt.Noisy        = t;
    tt.Filtered     = t;
    tt.LateNoise    = t;
    tt.NonLinearity = t;
    tt.DownSampled1 = downsample(t, downSampleFactor(1));
    tt.DownSampled2 = downsample(t, downSampleFactor(2));
    tt.DownSampled3 = downsample(t, downSampleFactor(3));
    
    x = rmfield(x, {'Template'});     
    
    datatypes = ['Ideal'; fieldnames(x)];
    
    subplot(221); cla;  hold on; plot(t, signal1, 'k-', 'LineWidth', 4); title('Signal A', 'FontSize', 20)
    subplot(222); cla;  hold on; plot(t, signal2, 'k-', 'LineWidth', 4); title('Signal B', 'FontSize', 20)

    for jj = 1:length(datatypes)
                
        if contains(datatypes{jj}, 'Ideal')
            LL1 = sum(log(normpdf(1/earlyNoiseLevel*(x.Noisy - signal1))));
            LL2 = sum(log(normpdf(1/earlyNoiseLevel*(x.Noisy - signal2))));
            Response = LL1 > LL2;
            PercentCorrect(jj,ii) = mean(Response' == labels)*100;
            
        else
            data = x.(datatypes{jj});
            cvmdl = fitcsvm(data', labels, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
            classLoss = kfoldLoss(cvmdl);
            PercentCorrect(jj,ii) = (1-classLoss) * 100;
            y1 = mean(data(:,find(labels==0,1)),2);
            y2 = mean(data(:,find(labels==1,1)),2);
            
            subplot(221), plot(tt.(datatypes{jj}), y1+jj*earlyNoiseLevel*2,'LineWidth', 2); ylim(earlyNoiseLevel*2*[-.5 length(datatypes)+1])
            subplot(222), plot(tt.(datatypes{jj}), y2+jj*earlyNoiseLevel*2,'LineWidth', 2); ylim(earlyNoiseLevel*2*[-.5 length(datatypes)+1])
            drawnow();
        end
    end
    
    legend(['noiseless'; datatypes(2:end)], 'Location', 'best'); 
        
   
    disp(PercentCorrect(:,ii)')
           
    subplot(2,2,3:4); cla
    bar(PercentCorrect'); set(gca, 'FontSize', 12)
    xlabel('Experiment number')
    ylabel('Accuracy')
    ylim([50 100])
    legend(datatypes, 'Location', 'eastoutside');


end


function dataout = applynonlinearity(datain, which_nonlin)
    
switch which_nonlin
    case 'ThresholdAndSquare'
        data  = abs(zscore(datain));
        thresh = .5;
        data(abs(data)<thresh) = 0;
        dataout = data.^2;

    case 'normalization'
        data = abs(zscore(datain));
        sigma = .3;
        n = 2;
        dataout = (data.^n)./(sigma.^n + data.^n);
        
    case 'none'
        dataout = datain;
end

end
