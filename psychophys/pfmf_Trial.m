function [trial, data] = pfmf_Trial(display, stimParams, data)
% trial = pfmf_Trial(display, stimParams)

numStimColors = display.stimRgbRange(2)-display.stimRgbRange(1)+1;

% We want to guarantee that the stimulus is modulated about the background
midStimColor = display.backColorIndex;

% display.gamma contains the gamma values needed to achieve each of
% numGammaEntries luminances, in a linear map. Thus, numGammaEntries/2 will
% yield the mean luminance, 1 will yield the minimum luminance, and
% numGammaEntries will yield the max luminance.
numGammaEntries = size(display.gamma,1);
midGammaIndex   = round(numGammaEntries/2);
halfStimRange   = stimParams.contrast*(numGammaEntries-1)*0.5;
gammaIndices    = linspace(-halfStimRange, halfStimRange, numStimColors)+midGammaIndex;

% if isfield(stimParams, 'present')
%     if ~stimParams.present, gammaIndices = gammaIndices * 0 + midGammaIndex; end
% end

cmap = zeros(display.maxRgbValue+1,3,1);

cmap(display.stimRgbRange(1)+1:display.stimRgbRange(2)+1,:) = display.gamma(round(gammaIndices),:);
% poke in reserved colors
cmap(1,:) = [0 0 0];
cmap(end,:) = [1 1 1];
% poke in the exact background color so that we can be sure it is in the
% center of the cmap (it might not be due to rounding).
cmap(display.backColorIndex, :) = display.gamma(midGammaIndex, :);

% Create the sequence
numFrames = round(display.frameRate * stimParams.duration);
seq = [-1 (1:numFrames+1)];

% Compute the images (if needed)
%   only compute images if they don't exist yet (we aren't manipulating
%   anything that requires recomputing the images, only the colormaps need
%   to be recomputed).
% if(~exist('data','var') | isempty(data))
% Generate images for both positions
radiusPix =  2*floor(angle2pix(display, stimParams.size/2)/2)+1;
spreadPix =  2*floor(angle2pix(display, stimParams.spread)/2)+1; % this will ensure that spreadPix is an odd number

[x,y] = meshgrid(-radiusPix:radiusPix,-radiusPix:radiusPix);
sz = size(x);

switch lower(stimParams.temporalEnvelopeShape)
    case 'gaussian'
        t = stimParams.temporalSpread/stimParams.duration;
        temporalWindow = exp(-.5*(([.5:numFrames]-numFrames/2)./(t*numFrames)).^2);
    case 'raisedcos'
        temporalWindow = ones(1,numFrames);
        len = ceil((stimParams.duration-stimParams.temporalSpread)/2*display.frameRate);
        endWin = (cos([0:len]/len*pi)+1)/2;
        temporalWindow(end-len+1:end) = endWin;
        temporalWindow(1:len) = fliplr(endWin);
    case 'flicker' % 1/2 cycle
        temporalWindow = sin(pi/numFrames*(1:numFrames));        
    otherwise
        temporalWindow = ones(1,numFrames);
end

sf       = stimParams.cyclesPerDegree*display.pixelSize*2*pi;
phaseInc = stimParams.cyclesPerSecond/display.frameRate*2*pi;
angle    = stimParams.orientDegrees*pi/180;
a        = cos(angle)*sf;
b        = sin(angle)*sf;
img      = cell(numFrames+1,1);
phase    = stimParams.phaseDegrees*pi/180;
spatialWindow = exp(-((x/spreadPix).^2)-((y/spreadPix).^2));

for ii=1:numFrames
    phase = phase+phaseInc;
    img{ii} = temporalWindow(ii)*spatialWindow.*sin(a*x+b*y+phase);
end
img{numFrames+1} = zeros(sz);
% compute grating
%img(:,:) = exp(-((x/spreadPix).^2) - ((y/spreadPix).^2)) ...
%			.* sin(x*.5*stimParams.cycles/radiusPix*2*pi);
% scale to the appropriate cmap range
for(ii=1:length(img))
    img{ii} = uint8(round(img{ii}.*(numStimColors/2-1)+midStimColor));
end

stimIm = createStimulusStruct(img, cmap, seq, []);
stimIm.imSize = sz;
clear('img');
stimIm = makeTextures(display, stimIm);
% else
% 	% the stimulus exists in data, so we just need to update the cmaps and seq
% 	data.cmap = cmap;
% 	data.s�eq = seq;
% end
% clear('cmap');
clear('seq');

sz = stimIm.imSize;
c = display.numPixels/2;
eccenPix = round(angle2pix(display, stimParams.eccentricity));
if stimParams.testPosition == 0 % West
    stimIm.destRect = round([c(1)-sz(1)/2 c(2)-sz(2)/2-eccenPix]);
elseif stimParams.testPosition == 90 % South
    stimIm.destRect = round([c(1)-sz(1)/2-eccenPix c(2)-sz(2)/2]);
elseif stimParams.testPosition == 180 % East
    stimIm.destRect = round([c(1)-sz(1)/2 c(2)-sz(2)/2+eccenPix]);  
elseif stimParams.testPosition == 270 % North
    stimIm.destRect = round([c(1)-sz(1)/2+eccenPix c(2)-sz(2)/2]);
else % fovea
    stimIm.destRect = round([c(1)-sz(1)/2 c(2)-sz(2)/2]);
end

stimIm.destRect = [stimIm.destRect stimIm.destRect+sz];

% get blank
blankIm = uint8(round(stimIm.imSize.*midStimColor));

% get preblank images: (blank with border)
numFramesPre = round(display.frameRate * stimParams.preStimduration);
seqPre = [-1 (1:numFramesPre+1)];
imgPre{numFramesPre+1} = zeros(sz);

blankCueIm = uint8(round(stimIm.imSize.*midStimColor));

for ii = 1:length(imgPre)
    
    % Create cue border
    blankCueIm([c(1)-sz(1)/2 c(2)-sz(2)/2-5]) = zeros(1,10);
    blankCueIm(borderLoc(3):borderLoc(2)+10,:)   = zeros(1,20);
    blankCueIm(borderLoc(3):borderLoc(3)-20,:)   = zeros(1,20);
    blankCueIm(borderLoc(4):borderLoc(4)-20,:)   = zeros(1,20);

    
    imgPre{ii} = blankCueIm;
end

preIm = createStimulusStruct(imgPre, cmap, seqPre, []);
preIm = makeTextures(display, preIm);

% Get sequence
if isfield(stimParams, 'present')
    % if stimulus is not present, make sequence a series of blank screens
    % by showing final (blank) screen in stimulus for every frame
    if ~stimParams.present, stimIm.seq = stimIm.seq * 0 + stimIm.seq(end); end
    stimIm.seq = stimIm.seq;
end

%% Make the trial
trial = addTrialEvent(display,[],'stimulusEvent', 'stimulus', preIm); %, 'duration', fixPeriodB4Stim);
trial = addTrialEvent(display,trial,'stimulusEvent', 'stimulus', stimIm);
trial = addTrialEvent(display,trial,'ITIEvent', 'stimulus', blankIm, 'duration', stimParams.iti);

% fixSeq   = ones(size(sequence));
% isi.sound = soundFreqSweep(500, 1000, .01);
% trial = addTrialEvent(display,[],'soundEvent', isi );
% isi.sound = zeros(1,1000);
% trial = addTrialEvent(display,trial,'soundEvent', isi );
% trial = addTrialEvent(display,trial,'stimulusEvent', 'stimulus', stimIm);

sequence      = [preIm.seq, stimIm.seq];
num_images    = length(imgPre)+length(imgStim);
timing        = (1:num_images) / display.frameRate ;
fixSeq        = ones(size(sequence));


data = 'done';

return;



% % specify the sequence of images as a vector of image indices
% sequence = [1 2 3]; % blank+cue stimulus blank
% timing   = (0:2) * 1/display.frameRate * stimParams.stimDuration; %stimParams.stimDuration;
% fixSeq   = ones(size(sequence));
% 
% fixPeriodB4Stim = 0.5; % second


end
