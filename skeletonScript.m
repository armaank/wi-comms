% A skeleton BER script for a wireless link simulation
clear all; close all; clc
numIter = 100; % The number of iterations of the monte-carlo simulation
nSym = 1000; % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = 2; % The M-ary number, 2 cooresponds to binary modulation

chan = 1; % No channel
%chan = [.1 .2 .4]; % Somewhat invertible channel, moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227'; % Not so invertible, severe ISI

% Time-varying Rayleigh multipath channel, try it if you dare. Or take
% wireless comms next semester.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to plot(chan)

% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

% Run the simulation numIter amount of times
for ii = 1:numIter
    
    bits = randi([0 log2(M)],1,nSym*log2(M)); % Generate random bits
    % New bits must be generated at every iteration
    
    % If you increase the M-ary number, you'll need to convert the its to
    % integers. See BIN2DE function
    % for binary, our MSG signal is simply the bits
    msg = bits;
    
    for jj = 1:lenSNR % one iteration of the simulation at each SNR value
        
        tx = qammod(msg, M); % BPSK modulate the signal
        
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx); % Apply the channel
        end
        
        txNoisy = awgn(txChan,SNR_Vec(jj),'measured'); % Add AWGN
        
        rx = qamdemod(txNoisy, M); % Demodulate
        
        % Again, if M was a larger number, need to convert symbols back to
        % bits here
        rxMSG = rx;
        
        % Compute and store the BER for this iteration
        
        [~, berVec(ii,jj)] = biterr(msg,rxMSG); % We're interested in the 
        % BER, which is the second output of BITERR
    
    end % End SNR iteration
end % End numIter iteration

% Compute and plot the mean BER
ber = mean(berVec, 1);

semilogy(SNR_Vec,ber)

% Compute the theoretical BER for this scenario, needs to be changed for
% different modulation schemes. There is no theoretical BER for multipath
% channels
berTheory = berawgn(SNR_Vec,'psk',2,'nondiff');
hold on;
semilogy(SNR_Vec+3,berTheory,'r')
legend('BER','Theoretical BER')








