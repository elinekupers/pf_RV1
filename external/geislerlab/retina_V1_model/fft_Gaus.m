function out = fft_Gaus(std_x, freq)

%fft_Gaus gives us the Fourier transform of a Gaussian with standard deviation std_x at a frequency freq.

std_f = 1/(2*pi*std_x);
out = exp(-(1/2)*(freq.^2)./(std_f.^2));