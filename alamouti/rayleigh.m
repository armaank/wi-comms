% simulating a rayleigh flat-fading channel
clc; clear all;
% generating fading spectrum

% params
fd = 10; % [hz] the doppler spread
fc = 0; % [hz] center frequency - might not need this ie fc = 0
N = 64; % number of samples in spectra

% generating conjugate-symetric random gaussian vectors
% number of points depends on the sampling frequency (spacing between
% points)
% white gaussian noise, 0dB variance 
v1 = 1/sqrt(2)*[randn(1,N/2) + j*randn(1,N/2)];
v1 = [fliplr(conj(v1)), v1];
v2 = 1/sqrt(2)*[randn(1,N/2) + j*randn(1,N/2)];
v2 = [fliplr(conj(v2)), v2];

% Computing doppler spectra
f = linspace(fc-fd, fc+fd, N)
% compute sampling time
ts = f(2) - f(1);
% Jakes doppler spectrum
S = (1.5)./(pi*fd*sqrt(1-((f-fc)/fd).^2))
% issue is that at bdry pts, S is inf
S(isinf(S)|isnan(S)) = 0; % Replace NaNs and infinite values with zeros

% this is the square of the desired frequency response of our doppler
% fading spectra
S = sqrt(S)

% since a gaussian random var in time is also a gaussian random var in
% frequency, we can simply multiply the Jakes spectra by our random
% variables and ifft to filter our time domain signal

v1_filtered = ifft(S.*v1, 2*N)
v2_filtered = ifft(S.*v2, 2*N)

chan = sqrt( (v1_filtered).^2 + (v2_filtered).^2)
figure(2)
plot(10*log10(abs(chan)))


