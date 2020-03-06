% simulating a rayleigh flat-fading channel

function chan = rayleigh(fd, N)

% generating conjugate-symetric random gaussian vectors
N = N/2;
v1 = sqrt(1/2)*[randn(1,N/2) + 1j*randn(1,N/2)];
v1 = [fliplr(conj(v1)), v1];
v2 = sqrt(1/2)*[randn(1,N/2) + 1j*randn(1,N/2)];
v2 = [fliplr(conj(v2)), v2];

% computing doppler spectra
f = linspace(-fd, fd, N);
S = (1.5)./(pi*fd*sqrt(1-((f)/fd).^2));
% issue is that at bdry pts, S is inf, sub NaNs and inf values with zero

S(isinf(S)|isnan(S)) = 0; 

S = sqrt(S);

v1_filtered = ifft(S.*v1, 2*N);
v2_filtered = ifft(S.*v2, 2*N);

chan = transpose(sqrt( (v1_filtered).^2 + (v2_filtered).^2));

end
