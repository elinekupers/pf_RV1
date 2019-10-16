function out = MTF_optics(freq)

%MTF_optics codes the modulation transfer function (MTF) for the optics of the human eye. The input can be an entire
%matrix of frequencies instead of just a single frequency. Equation and parameters are from Navarro et al. (1993).

out = (1-0.22).*exp(-0.172.*freq) + 0.22.*exp(-0.037.*freq); %for the fovea



