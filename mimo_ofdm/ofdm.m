%% OFDM simulation
clc; clear all;
%% params

nIter = 3;  % The number of iterations of the simulation

N = 64; % number of subcarriers
mu = 16; % length of cyclic prefix
nSamp = N + mu; % total number of samples per ofdm symbol period
M = 16; % modulation order for QAM
npt = 64; % number of fft/ifft points


chan1 = [.1 .2 .4]; % Somewhat invertible channel, moderate ISI
chan2 = [0.227 0.460 0.688 0.460 0.227]; % Not so invertible, severe ISI

EbNo = -5:2:20; % 
snr_vec = EbNo + 10*log10(N/nSamp); % snr conversion from EbNo
lenSNR = length(snr_vec); 

ber = zeros(2,lenSNR,nIter);

%% wls - ofdm zf 

for ii = 1:nIter
    
%     bits = randi([0,M-1],1, nSym*(N-mu));     % Generate random bits
    mod_data = qammod(bits,M);  % modulate the signal
    
end

%% plot

ber = mean(ber,nIter); %take mean across all iterations

%plotting ber for different equalization techniques
semilogy(EbNo,ber(1,:),'DisplayName','OFDM zf');
hold on;

title('BER curves');
xlabel('EbNo(dB)');
ylabel('Bit Error Rate');
legend('show');






