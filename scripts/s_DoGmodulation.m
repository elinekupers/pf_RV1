%% Effect of DoG filters on Gabor stimulus

% Mosaic params
ncones = 79;
fov = 2; % degrees
x = (0.5:ncones-.5)/ncones*fov; % degrees

% RGC params kaplan
szratio = 6;
wc = .64;
ws = 1-wc;

% RGC params bradley abrams geisler 2014
% szratio = 10;
% wc = 0.53;
% ws = 1-wc;

saveDir = fullfile(pfRV1rootPath, 'figures', 'DoG1D');
if ~exist('saveDir', 'dir'); mkdir(saveDir); end

%% Load or simulate stimulus
load(fullfile(pfRV1rootPath, 'external', 'data', 'stimulus_vs_noise1'), 'noisyImage', 'theNoise', 'absorptions')
stim2D      = mean(absorptions,3);
noise2D     = mean(theNoise,3);
coneresp2D  = mean(noisyImage,3);
midpoint    = ceil(size(stim2D,1)/2);

figure(100); set(gcf, 'Position', [16 487 1641 461]); clf;
subplot(131); imagesc(stim2D); colormap gray; 
hold on; plot([1 size(stim2D,1)],[midpoint, midpoint],  'r-', 'LineWidth',3); set(gca, 'TickDir', 'out', 'FontSize', 14)
xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off;
title('Stimulus (L-cone only - contrast: 1% - sf: 4 cpd - fov: 2 deg')

subplot(132); imagesc(noise2D); colormap gray; set(gca, 'TickDir', 'out', 'FontSize', 14)
hold on; plot([1 size(stim2D,1)],[midpoint, midpoint],  'r-', 'LineWidth',3);
xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off;
title('Noise')

subplot(133); imagesc(coneresp2D); colormap gray; set(gca, 'TickDir', 'out', 'FontSize', 14)
hold on; plot( [1 size(stim2D,1)],[midpoint, midpoint], 'r-', 'LineWidth',3);
xlabel('Cone (columns)'), xlabel('Cone (rows)'); box off;
title('Cone absorptions (stim+noise)')

print(gcf, fullfile(saveDir, 'Stim_Noise_ConeResp'), '-dpng')


% Go 1D for now
stim        = stim2D(midpoint,:)';
noise       = noise2D(midpoint,:)';
coneresp    = coneresp2D(midpoint,:)';

stimfreq      = 4; % cpd
sigmastim     = .5; % degrees
%
% % build stim
% stimgauss     = exp(-(x-fov/2).^2 ./ (2*sigmastim^2));
% stimgauss     = 2*stimgauss / sum(stimgauss);
% stimgrating   = sin(stimfreq * x * 2 * pi);
% stim          = stimgauss .* stimgrating;
% noise         = randn(ncones,1) * 0.01 * ones(1,5);
% coneresp  = noise + stim';

% Stim and noise in FFT
STIM     = abs(fft(stim(:)));
CONERESP = abs(fft(coneresp(:)));
NOISE    = abs(fft(noise(:)));

hc = NaN(size(coneresp));
hs = NaN(size(coneresp));

numC2RGCratios = 5;

cmap = lines(numC2RGCratios);

for ii = 1:numC2RGCratios
    
    conespacing = 1/ncones*fov;
    
    sigmac = ii/3 * conespacing;
    
    sigmas = sigmac * szratio;
    
    
    % Get DoG center, and surround
    gc = 1/(sqrt(2*pi)*sigmac) * exp(-(x-fov/2).^2 ./ (2*sigmac^2));
    
    gs = 1/(sqrt(2*pi)*sigmas) * exp(-(x-fov/2).^2 ./ (2*(sigmas).^2));
    
    gc = gc/sum(gc);
    
    gs = gs/sum(gs);
    
    hc(:,ii) = 2*wc*gc;
    hs(:,ii) = 2*ws*gs;
    labels{ii} = sprintf('Cone:RGC=1:%d',ii);
end

f  = hc - hs;

rgcc = conv2(coneresp, hc, 'same');
rgcs = conv2(coneresp, hs, 'same');
rgcresp = rgcc - rgcs; % conv(coneresp, f, 'same');


Gc = abs(fft(hc));
Gs = abs(fft(hs));
G  = abs(fft(f));

CONESOUTPUT = bsxfun(@times, G, CONERESP);
STIMOUTPUT  = bsxfun(@times, G, STIM);
NOISEOUTPUT = bsxfun(@times, G, NOISE);

fs = (0:ncones-1)/2;

figure(1), clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);


subplot(4,2,1)
plot(x, hc, 'LineWidth', 4); 
ylim([0 1.25]); title ('Center'); 
box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('modulation (a.u.)'); 

subplot(4,2,3)
plot(x, hs, 'LineWidth', 4); 
ylim([0 1.25]); title ('Surround'); 
box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('modulation (a.u.)'); 
legend(labels); legend boxoff;

subplot(4,2,5)
plot(x, f, 'LineWidth', 4); 
ylim([-0.2 1.25]); title ('DoG');
box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('modulation (a.u.)'); 

subplot(4,2,7)
plot(x, stim-min(stim), 'k-', x, noise-min(noise), 'k--', 'LineWidth',2); 
set(gca, 'TickDir', 'out', 'FontSize', 14)
title ('Stim + noise (scaled)'); 
box off; legend({'Stimulus', 'Noise'}); legend boxoff;
ylabel('absorption (count)'); xlabel('Visual field (deg)')

subplot(4,2,2)
plot(fs, Gc, 'LineWidth', 4); 
xlim([0 max(fs)/2]); ylim([0 1.5]);
title ('Center'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('Amplitude')

subplot(4,2,4)
plot(fs, Gs, 'LineWidth', 4);
xlim([0 max(fs)/2]); ylim([0 1.5]);
title ('Surround'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('Amplitude')
legend(labels); legend boxoff;

subplot(4,2,6)
plot(fs, G, 'LineWidth', 4); 
xlim([0 max(fs)/2]); ylim([0 1.5]); 
title ('DoG'); box off; set(gca, 'TickDir', 'out', 'FontSize', 14)
ylabel('Amplitude')

subplot(4,2,8)
plot(fs, STIM-min(STIM), 'k-', fs, NOISE-min(NOISE), 'k--', 'LineWidth', 2);
xlim([0 max(fs)/2]); ylim([0 50]); 
title ('Stim + Noise (scaled)'); 
legend({'Stimulus', 'Noise'}); legend boxoff; 
set(gca, 'TickDir', 'out', 'FontSize', 14)
xlabel('Spatial frequency (cpd)'); ylabel('Amplitude'); box off;

print(gcf, fullfile(saveDir, 'DoG1D_filter_stim'), '-dpng')



figure(2);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
for ii = 1:numC2RGCratios
    
    % Stim (no noise)
    subplot(311); title('FFT RGC outputs for STIM');
    plot(fs, STIMOUTPUT(:,ii), '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    xlim([0 max(fs)/2]); ylim([0 50]); set(gca, 'TickDir', 'out', 'FontSize', 14);
    legend(labels); legend boxoff; ylabel('Amplitude')
    box off;
    
    % Noise
    subplot(312); title('FFT RGC outputs for NOISE');
    plot(fs, NOISEOUTPUT(:,ii), '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    xlim([0 max(fs)/2]); ylim([0 50]); set(gca, 'TickDir', 'out', 'FontSize', 14)
    box off; ylabel('Amplitude')
    
    % Cones (stim+noise)
    subplot(313); title('FFT RGC outputs for CONES (stim+noise)');
    xlim([0 max(fs)/2]); ylim([0 50]);
    plot(fs, CONESOUTPUT(:,ii), '-','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    xlabel('Spatial frequency (cpd)'); set(gca, 'TickDir', 'out', 'FontSize', 14);
    box off; ylabel('Amplitude')
end

print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_individual'), '-dpng')


figure(3);  clf; set(gcf, 'Color', 'w', 'Position', [1 117 1056 831]);
for ii = 1:numC2RGCratios
    
    % Stim (no noise)
    title('FFT RGC outputs');
    plot(fs, STIMOUTPUT(:,ii), '-', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    plot(fs, NOISEOUTPUT(:,ii), '--', 'LineWidth', 2, 'Color', cmap(ii,:)); hold on;
    plot(fs, CONESOUTPUT(:,ii), ':','LineWidth', 2, 'Color', cmap(ii,:)); hold on;
end

xlim([0 max(fs)/2]); ylim([0 50]); set(gca, 'TickDir', 'out', 'FontSize', 14);
legend({'Stim', 'Noise', 'Cones'}); legend boxoff;
ylabel('Amplitude')
xlabel('Spatial frequency (cpd)');
box off;

print(gcf, fullfile(saveDir, 'DoG1D_RGCOUTPUT_all'), '-dpng')

return
