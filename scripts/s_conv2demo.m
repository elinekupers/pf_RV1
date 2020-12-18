%% Effect of zero-padding using conv2 to filter mRGC responses

% general flags
useArtificalStim = true;
doSubsampling    = false;
doSmallerArray   = false;

% Mosaic params
ncones    = 79;
conearray = zeros(ncones, ncones);
fov       = 2; % degrees
x = (0.5:ncones-.5)/ncones*fov; % degrees

% RGC params kaplan
szratio = 6;
wc = .64;
ws = 1-wc;

saveFig = false;
if saveFig
    saveDir = fullfile(pfRV1rootPath, 'figures', 'conv2demo1');
    if ~exist('saveDir', 'dir'); mkdir(saveDir); end
end
%% Load or simulate stimulus

% Stim params
stimfreq      = 4; % cpd
sigmastim     = .5; % degrees

% Build stim if requested
if useArtificalStim
    stimgauss     = exp(-(x-fov/2).^2 ./ (2*sigmastim^2));
    stimgauss     = 2*stimgauss / sum(stimgauss);
    stimgrating   = sin(stimfreq * x * 2 * pi);
    stim          = stimgauss .* stimgrating;
    stim          = stim';
    noise         = randn(ncones,1) * 0.01 * ones(1,5);
    coneresp      = noise + stim;
    
     % Get snr
    snr_stim     = snr(stim-mean(stim), noise(:,1));
    snr_noise    = snr(noise-mean(noise), noise);
    snr_coneresp = snr(coneresp-mean(coneresp), noise);
    
else
    load(fullfile(pfRV1rootPath, 'external', 'data', 'stimulus_vs_noise1'), 'noisyImage', 'theNoise', 'absorptions')
    % Absorptions without noise (i.e. noiseless stimulus)
    stim2D      = mean(absorptions,3);
    % Absorption noise (i.e. no stimulus)
    noise2D     = mean(theNoise,3);
    % Absorptions with noise (i.e. stimulus + noise)
    coneresp2D  = mean(noisyImage,3);
    
    % Get snr
    snr_stim     = snr(stim2D-mean(stim2D(:)), noise2D);
    snr_noise    = snr(noise2D-mean(noise2D(:)), noise2D);
    snr_coneresp = snr(coneresp2D-mean(coneresp2D(:)), noise2D);
    
    % Go 1D for now
    midpoint    = ceil(size(stim2D,1)/2);
    stim        = stim2D(midpoint,:)';% - mean(stim2D(midpoint,:));
    noise       = noise2D(midpoint,:)';
    coneresp    = coneresp2D(midpoint,:)';% - mean(coneresp2D(midpoint,:));
    
    % Plot it!
    figure(100); set(gcf, 'Position', [16 487 1641 461]); clf;
    subplot(131); imagesc(stim2D); colormap gray; axis square;
    hold on; plot([1 size(stim2D,1)],[midpoint, midpoint],  'r-', 'LineWidth',3); set(gca, 'TickDir', 'out', 'FontSize', 14)
    xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off; colorbar;
    title('Stimulus (L-cone only - contrast: 1% - sf: 4 cpd - fov: 2 deg')
    
    subplot(132); imagesc(noise2D); colormap gray; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
    hold on; plot([1 size(stim2D,1)],[midpoint, midpoint],  'r-', 'LineWidth',3);
    xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off; colorbar;
    title('Noise')
    
    subplot(133); imagesc(coneresp2D); colormap gray; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
    hold on; plot( [1 size(stim2D,1)],[midpoint, midpoint], 'r-', 'LineWidth',3);
    xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off; colorbar;
    title('Cone absorptions (stim+noise)')
    
    if saveFig
        print(gcf, fullfile(saveDir, 'Stim_Noise_ConeResp'), '-dpng')
    end
    
end

%% Make DoG filters

hc = NaN(size(coneresp));
hs = NaN(size(coneresp));

numC2RGCratios = 5;

cmap = lines(numC2RGCratios);

for ii = 1:numC2RGCratios
    
    % Create resampling vector with cone2RGC ratio.
    rowIndices{ii} = 1:ii:ncones;
    
    % Get cone spacing
    conespacing = 1/ncones*fov;
    
    % Set RGC center and surround sizes
    sigmac = ii/3 * conespacing;
    sigmas = sigmac * szratio;
    
    % Get DoG center, and surround
    gc = 1/(sqrt(2*pi)*sigmac) * exp(-(x-fov/2).^2 ./ (2*sigmac^2));
    gs = 1/(sqrt(2*pi)*sigmas) * exp(-(x-fov/2).^2 ./ (2*(sigmas).^2));
    
    % Normalize area
    gc = gc/sum(gc);
    gs = gs/sum(gs);
    
    % Weight center/surround
    hc(:,ii) = 2*wc*gc;
    hs(:,ii) = 2*ws*gs;
    
    % Get labels for plots
    labels{ii} = sprintf('RGC:Cones=1:%d',ii);
    
    if doSmallerArray
        % get reasonable estimate how many cones the surround Gaussian captures,
        % skip these cones on either side to avoid summing over zero padded
        % part
        startArray = ceil(max(diff(find(gs<0.0001)))/2);
        subtractAtEndArray = startArray;
        rowIndiced_smallerArray{ii} = startArray:ii:(ncones-subtractAtEndArray);
    end
end

% Get DoG filter by subtracting center and surround gaussian
f  = hc - hs;

%% Filter the cone absorptions to get mRGC responses

coneresp = repmat(coneresp, [1,5]);
stim     = repmat(stim, [1,5]);
noise    = repmat(noise, [1,5]);

rgcresp = [];
rgcstim = rgcresp;
rgcnoise = rgcresp;

for ii = 1:numC2RGCratios
    
    padMean1         = ones(ceil(length(coneresp(:,ii))/2)-1,1)*mean(coneresp(:,ii));
    coneresp_padMean(:,ii) = [iePoisson(padMean1); coneresp(:,ii); iePoisson(padMean1)];
    padMean2         = ones(ceil(length(stim(:,ii))/2)-1,1)*mean(stim(:,ii));
    stim_padMean(:,ii)     = [iePoisson(padMean2); stim(:,ii); iePoisson(padMean2)];
    padMean3         = ones(ceil(length(noise(:,ii))/2)-1,1)*mean(noise(:,ii));
    noise_padMean(:,ii)    = [iePoisson(padMean3); noise(:,ii); iePoisson(padMean3)];
        
    % Convolve cone response with DoG filter
    rgcresp(:,ii)  = conv(coneresp(:,ii), f(:,ii), 'same');
    rgcstim(:,ii)  = conv(stim(:,ii), f(:,ii), 'same');
    rgcnoise(:,ii) = conv(noise(:,ii), f(:,ii), 'same');
    
    % Convolve cone response with DoG filter
    rgcresp_padded_tmp  = conv(coneresp_padMean(:,ii), f(:,ii), 'same');
    rgcstim_padded_tmp  = conv(stim_padMean(:,ii), f(:,ii), 'same');
    rgcnoise_padded_tmp = conv(noise_padMean(:,ii), f(:,ii), 'same');
    
    rgcresp_padded(:,ii)  = rgcresp_padded_tmp(length(padMean1)+1:end-length(padMean1));
    rgcstim_padded(:,ii)  = rgcstim_padded_tmp(length(padMean2)+1:end-length(padMean2));
    rgcnoise_padded(:,ii) = rgcnoise_padded_tmp(length(padMean3)+1:end-length(padMean3));
    
    
    snrFilteredAbsorptions(:,ii) = snr(rgcresp(:,ii)-mean(rgcresp(:,ii)), rgcnoise(:,ii));
    snrFilteredAbsorptions_padded(:,ii) = snr(rgcresp_padded(:,ii)-mean(rgcresp_padded(:,ii)), rgcnoise_padded(:,ii));
    
    
    if doSubsampling
        
        if doSmallerArray
            % No padding: Smaller array for resampling
            rgcresp_sub{ii}  = rgcresp(rowIndiced_smallerArray{ii},ii);
            rgcstim_sub{ii}  = rgcstim(rowIndiced_smallerArray{ii},ii);
            rgcnoise_sub{ii} = rgcnoise(rowIndiced_smallerArray{ii},ii);
            
            RGCOUTPUT{ii} = abs(fft(rgcresp_sub{ii}));
            STIMOUTPUT{ii}  = abs(fft(rgcstim_sub{ii}));
            NOISEOUTPUT{ii} = abs(fft(rgcnoise_sub{ii}));
            
            % With padding: Smaller array for resampling
            rgcresp_sub_padded{ii}  = rgcresp(rowIndiced_smallerArray{ii},ii);
            rgcstim_sub_padded{ii}  = rgcstim(rowIndiced_smallerArray{ii},ii);
            rgcnoise_sub_padded{ii} = rgcnoise(rowIndiced_smallerArray{ii},ii);
            
            RGCOUTPUT_padded{ii} = abs(fft(rgcresp_sub_padded{ii}));
            STIMOUTPUT_padded{ii}  = abs(fft(rgcstim_sub_padded{ii}));
            NOISEOUTPUT_padded{ii} = abs(fft(rgcnoise_sub_padded{ii}));
        else      
            %Subsample responses without padding
            rgcresp_sub{ii}  = rgcresp(rowIndices{ii},ii);
            rgcstim_sub{ii}  = rgcstim(rowIndices{ii},ii);
            rgcnoise_sub{ii} = rgcnoise(rowIndices{ii},ii);
            
            % FFT: WITH subsampling (after filter)
            RGCOUTPUT{ii} = abs(fft(rgcresp_sub{ii}));
            STIMOUTPUT{ii}  = abs(fft(rgcstim_sub{ii}));
            NOISEOUTPUT{ii} = abs(fft(rgcnoise_sub{ii}));
            
            % Subsample responses with padding
            rgcresp_sub_padded{ii}  = rgcresp_padded(rowIndices{ii},ii);
            rgcstim_sub_padded{ii}  = rgcstim_padded(rowIndices{ii},ii);
            rgcnoise_sub_padded{ii} = rgcnoise_padded(rowIndices{ii},ii);
            
            % and their FFT
            RGCOUTPUT_padded{ii} = abs(fft(rgcresp_sub_padded{ii}));
            STIMOUTPUT_padded{ii}  = abs(fft(rgcstim_sub_padded{ii}));
            NOISEOUTPUT_padded{ii} = abs(fft(rgcnoise_sub_padded{ii}));
        end
        
    else
        % FFT: NO subsampling (just filter)
        RGCOUTPUT{ii} = abs(fft(rgcresp(:,ii)));
        STIMOUTPUT{ii}  = abs(fft(rgcstim(:,ii)));
        NOISEOUTPUT{ii} = abs(fft(rgcnoise(:,ii)));
        
        RGCOUTPUT_padded{ii} = abs(fft(rgcresp_padded(:,ii)));
        STIMOUTPUT_padded{ii}  = abs(fft(rgcstim_padded(:,ii)));
        NOISEOUTPUT_padded{ii} = abs(fft(rgcnoise_padded(:,ii)));
        
    end
end
    % Stim and noise in FFT
    STIM     = abs(fft(stim));
    CONERESP = abs(fft(coneresp));
    NOISE    = abs(fft(noise));
    
    % Get FFT amplitudes of Gaussian filters
    Gc = abs(fft(hc));
    Gs = abs(fft(hs));
    G  = abs(fft(f));
    
    % % Multiply with FFT of cone/stim/noise responses
    % CONESOUTPUT = bsxfun(@times, G, CONERESP);
    % STIMOUTPUT  = bsxfun(@times, G, STIM);
    % NOISEOUTPUT = bsxfun(@times, G, NOISE);

    % FFT ampl x-axis
    fs = (0:ncones-1)/2;
    
    
    %% Visualize filters in visual and fft domain
    figure(1), clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    
    % Plot center Gaussian
    subplot(4,2,1)
    plot(x, hc, 'LineWidth', 4);
    ylim([0 1.25]); title ('Center');
    box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('modulation (a.u.)');
    
    % Plot surround Gaussian
    subplot(4,2,3)
    plot(x, hs, 'LineWidth', 4);
    ylim([0 1.25]); title ('Surround');
    box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('modulation (a.u.)');
    legend(labels); legend boxoff;
    
    % Plot surround DoG
    subplot(4,2,5)
    plot(x, f, 'LineWidth', 4);
    ylim([-0.2 1.25]); title ('DoG');
    box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('modulation (a.u.)');
    
    % Plot cone absorptions for stimulus + noise
    subplot(4,2,7)
    plot(x, stim, 'k-', x, noise, 'k--', 'LineWidth',2);
    set(gca, 'TickDir', 'out', 'FontSize', 14)
    title ('Stim + noise');
    box off; legend({'Stimulus', 'Noise'}); legend boxoff;
    ylabel('absorption (count)'); xlabel('Visual field (deg)')
    
    % Plot FFT center Gaussian
    subplot(4,2,2)
    plot(fs, Gc, 'LineWidth', 4);
    xlim([0 max(fs)/2]); ylim([0 1.5]);
    title ('Center'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('Amplitude')
    
    % Plot FFT surround Gaussian
    subplot(4,2,4)
    plot(fs, Gs, 'LineWidth', 4);
    xlim([0 max(fs)/2]); ylim([0 1.5]);
    title ('Surround'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('Amplitude')
    legend(labels); legend boxoff;
    
    % Plot FFT DoG
    subplot(4,2,6)
    plot(fs, G, 'LineWidth', 4);
    xlim([0 max(fs)/2]); ylim([0 1.5]);
    title ('DoG'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
    ylabel('Amplitude')
    
    % Plot FFT Stim and noise
    subplot(4,2,8)
    plot(fs, STIM, 'k-', fs, NOISE, 'k--', 'LineWidth', 2);
    xlim([0 max(fs)/2]); ylim([0 50]);
    title ('Stim + Noise');
    legend({'Stimulus', 'Noise'}); legend boxoff;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
    xlabel('Spatial frequency (cpd)'); ylabel('Amplitude'); box off;
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_filter_stim'), '-dpng')
    end
    %% Visualize inputs (noise, stim, noise+stim) as snr
    figure(12), clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    
    subplot(3,1,1)
    plot(x, stim, 'k', 'LineWidth',2);
    set(gca, 'TickDir', 'out', 'FontSize', 14); box off;
    title(sprintf('Stimulus absorptions SNR: %2.2f', snr_stim));
    ylabel('absorption (count)'); xlabel('Visual field (deg)')
    
    subplot(3,1,2)
    plot(x, noise, 'k', 'LineWidth',2);
    set(gca, 'TickDir', 'out', 'FontSize', 14); box off;
    title(sprintf('Noise absorptions SNR: %2.2f', snr_noise));
    ylabel('absorption (count)'); xlabel('Visual field (deg)')
    
    subplot(3,1,3)
    plot(x, coneresp, 'k', 'LineWidth',2);
    set(gca, 'TickDir', 'out', 'FontSize', 14); box off;
    title(sprintf('Total absorptions SNR: %2.2f', snr_coneresp));
    ylabel('absorption (count)'); xlabel('Visual field (deg)')
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_SNR_stim_noise'), '-dpng')
    end
    %%
    
    figure(13);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    for ii = 1:numC2RGCratios
        fs = (0:length(STIMOUTPUT{ii})-1)/2;
        yl = max(NOISEOUTPUT{ii})*5;
        % Stim (no noise)
        subplot(311); title('FFT RGC outputs for STIM');
        plot(fs, STIMOUTPUT{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14);
        legend(labels); legend boxoff; ylabel('Amplitude')
        box off;
        if ii==1; xlim([0 max(fs)/2]); end
        
        % Noise
        subplot(312); title('FFT RGC outputs for NOISE');
        plot(fs, NOISEOUTPUT{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14)
        box off; ylabel('Amplitude')
        if ii==1; xlim([0 max(fs)/2]); end
        
        % Cones (stim+noise)
        subplot(313); title('FFT RGC outputs for CONES (stim+noise)');
        ylim([0 yl]);
        plot(fs, RGCOUTPUT{ii}, '-','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        xlabel('Spatial frequency (cpd)'); set(gca, 'TickDir', 'out', 'FontSize', 14);
        box off; ylabel('Amplitude');
        if ii==1; xlim([0 max(fs)/2]); end
    end
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_fft'), '-dpng')
    end
    
    %%
    figure(14);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    for ii = 1:numC2RGCratios
        fs = (0:length(STIMOUTPUT_padded{ii})-1)/2;
        yl = max(NOISEOUTPUT_padded{ii})*5;
        % Stim (no noise)
        subplot(311); title('FFT RGC outputs for STIM,  w/ mean padding');
        plot(fs, STIMOUTPUT_padded{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14);
        legend(labels); legend boxoff; ylabel('Amplitude')
        box off;
        if ii==1; xlim([0 max(fs)/2]); end
        
        % Noise
        subplot(312); title('FFT RGC outputs for NOISE,  w/ mean padding');
        plot(fs, NOISEOUTPUT_padded{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14)
        box off; ylabel('Amplitude')
        if ii==1; xlim([0 max(fs)/2]); end
        
        % Cones (stim+noise)
        subplot(313); title('FFT RGC outputs for CONES (stim+noise),  w/ mean padding');
        ylim([0 yl]);
        plot(fs, RGCOUTPUT_padded{ii}, '-','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        xlabel('Spatial frequency (cpd)'); set(gca, 'TickDir', 'out', 'FontSize', 14);
        box off; ylabel('Amplitude');
        if ii==1; xlim([0 max(fs)/2]); end
    end
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_fft_meanpadding'), '-dpng')
    end
    
    %%
    figure(15);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    for ii = 1:numC2RGCratios
        fs = (0:length(STIMOUTPUT{ii})-1)/2;
        
        % Stim (no noise)
        title('FFT RGC outputs');
        plot(fs, STIMOUTPUT{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        plot(fs, NOISEOUTPUT{ii}, '--', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        plot(fs, RGCOUTPUT{ii}, ':','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    end
    yl = max(NOISEOUTPUT{ii})*20;
    xlim([0 18]); ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14);
    legend({'Stim', 'Noise', 'Cones'}); legend boxoff;
    ylabel('Amplitude')
    xlabel('Spatial frequency (cpd)');
    box off;
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_all'), '-dpng')
    end
    %%
    figure(16);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
    for ii = 1:numC2RGCratios
        
        fs = (0:length(STIMOUTPUT_padded{ii})-1)/2;
        
        % Plot FFT RGC
        title('FFT RGC outputs  w/ mean padding');
        plot(fs, STIMOUTPUT_padded{ii}, '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        plot(fs, NOISEOUTPUT_padded{ii}, '--', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
        plot(fs, RGCOUTPUT_padded{ii}, ':','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    end
     yl = max(NOISEOUTPUT{ii})*8;
    xlim([0 18]); ylim([0 yl]); set(gca, 'TickDir', 'out', 'FontSize', 14);
    legend({'Stim', 'Noise', 'Cones'}); legend boxoff;
    ylabel('Amplitude')
    xlabel('Spatial frequency (cpd)');
    box off;
    
    if saveFig
        print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_meanpadding'), '-dpng')
    end
    
    %%
    figure;
    for ii = 1:5
        subplot(5,1,ii);
        plot(coneresp(:,ii), 'lineWidth',2); hold all; plot(rgcresp(:,ii), 'lineWidth',2); plot(rgcresp_padded(:,ii), ':', 'lineWidth',2)
        set(gca, 'TickDir','out', 'FontSize',15)
        title(sprintf('Ratio %1.1f:1', 2/ii))
        
    end
    
    xlabel('cones')
    ylabel('absorptions')
    legend('cone absorptions', 'rgc response', 'rgc response w/ mean padding')
    
    if doSmallerArray
        figure;
        for ii = 1:5
            subplot(5,1,ii);
            x = 0:(size(coneresp,1)-1);
            x2 = 0:length(rgcresp_sub{ii})-1; 
            x2 = x2-(max(x2)/2);
            x2 = x2 + (length(x)/2);
            plot(x,coneresp(:,ii), 'lineWidth',2); hold all; 
            plot(x2, rgcresp_sub{ii}, 'lineWidth',2); plot(x2, rgcresp_sub_padded{ii}, ':', 'lineWidth',2)
            set(gca, 'TickDir','out', 'FontSize',15)
            title(sprintf('Ratio %1.1f:1', 2/ii))

        end

        xlabel('cones')
        ylabel('absorptions')
        legend('cone absorptions', 'rgc response', 'rgc response w/ mean padding')
    end
    
    
    return
