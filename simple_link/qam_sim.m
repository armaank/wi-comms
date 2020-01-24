%% ECE408 Wireless Communications
% Armaan Kohli
% Spring 2020
% Project 1
clear all;
close all;
clc;
%% QAM Simulation w/ No Channel

numIter = 5; % The number of iterations of the monte-carlo simulation
nSym = 1000; % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = [4, 16]; % The M-ary number

chan = 1; % No channel

% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

% Run the simulation numIter amount of times for both 4-ary and 16-ary QAM
for m = M
    for ii = 1:numIter
        
        K = log2(m); % COMMENT HERE
        % Generating random bits to send
        bits = randi(2,[nSym*K, 1])-1; 
        
        % Reshape bits into K-tuples
        data = reshape(bits,length(bits)/K,K);   
        msg = bi2de(data);

        for jj = 1:lenSNR % One iteration of the simulation at each SNR

            tx = qammod(msg,m); % QAM modulate the signal

            txChan = tx; % Perfect channel - introduces no ISI

            % Convert from EbNo to SNR and add AWGN
            txNoisy = awgn(txChan,(10*log10(K) + SNR_Vec(jj)),'measured'); 
            
            rx = qamdemod(txNoisy,m); % Demodulate

            % Converting symbols back to bits
            dataOutMatrix = de2bi(rx,K);
            rxMSG = dataOutMatrix(:); 

            % Compute and store the BER for this iteration

            [~, berVec(ii,jj)] = biterr(bits, rxMSG);  
        end % End SNR iteration
    end % End numIter iteration
    
    % Generate plots for each M-ary scheme
    
    % Generate scatterplot
    title1 = sprintf("Scatter Plot for %d - ary QAM with no ISI",m);
    sPlotFig = scatterplot(txNoisy,1,0,'g.');
    hold on;
    figure(1);
    scatterplot(tx,1,0,'k*',sPlotFig)
    hold off;
    title(title1)

    % Compute and plot the mean BER
    ber = mean(berVec,1);
    figure(2)
    semilogy(SNR_Vec,ber)
    % Compute the theoretical BER for this scenario and graph waterfall 
    berTheory = berawgn(SNR_Vec, 'qam', m); % theoretical BER
    title2 = sprintf("BER Waterfall Curve for %d - ary QAM with no ISI",m);
    hold on
    semilogy(SNR_Vec, berTheory)
    xlabel('SNR(dB)');  ylabel('BER');
    legend('BER','Theoretical BER','Location','southwest')
    title(title2)

end % End M-ary iteration

%% Part 2: Equalizing a Channel with Moderate ISI for BPSK

numIter = 1; 
nSym = 100000;    
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = 2;        % The M-ary number, 2 corresponds to binary modulation
K = log2(M);  % Number of bits per symbol
trainlength = nSym/10;

chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI

% Modulator
bpskMod = comm.BPSKModulator('PhaseOffset',0);
% Demodulator
bpskDemod = comm.BPSKDemodulator('PhaseOffset',0); 

% Equalizer Object
eq1 = dfe(2,3,signlms(0.01));
eq1.SigConst = bpskMod((0:M-1)')'; % Set signal constellation.
eq1.ResetBeforeFiltering=0;

berVec1 = zeros(numIter, lenSNR);

for i = 1:numIter
   
    bits = randi(2,[nSym*K, 1])-1;
    data = reshape(bits,length(bits)/K,K);   % Reshape bits into binary k-tuples, K = log2(M)

    msg = bi2de(data);
    
    for j = 1:lenSNR 
        
        tx = bpskMod(msg);
        
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
        
        % Convert from EbNo to SNR.
        % Note: Because No = 2*noiseVariance^2, we must add ~3 dB
        % to get SNR (because 10*log10(2) ~= 3).
        txNoisy = awgn(txChan,(10*log10(K)+SNR_Vec(j)),'measured'); % Add AWGN
        
        EqData1 = equalize(eq1,txNoisy,tx(1:trainlength));
        rx1 = bpskDemod(EqData1);
                             
        dataOutMatrix1 = de2bi(rx1,K);
        rxMSG1 = dataOutMatrix1(:); 
               
        [~, berVec1(i,j)] = biterr(bits(2:end), rxMSG1(2:end));  % We're interested in the BER, which is the 2nd output of BITERR
    end
    
end

h = scatterplot(txNoisy,1,trainlength,'bx'); hold on;
scatterplot(EqData1,1,trainlength,'g.',h);
scatterplot(eq1.SigConst,1,0,'k*',h);
legend('Uneqalized signal','Equalized signal',...
   'Ideal signal constellation');
hold off;


ber1 = mean(berVec1,1);
figure(2)
semilogy(SNR_Vec, ber1)
hold on;
berTheory = berawgn(SNR_Vec,'psk',M,'nondiff');
semilogy(SNR_Vec,berTheory,'r')
%line(12,,'Color','red','LineStyle','--');
legend('BER','Theoretical BER')








