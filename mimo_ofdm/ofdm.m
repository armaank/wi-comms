%% OFDM simulation
clc; clear all;
%% simulation parameters

% modulation order
M = 16;
K = log2(M);     

% symbols per transmission 
N_syms = 1e4;

% number of channels
n_chan = 3;

EbNo = -10:2:30;                  
snr = EbNo + 10*log10(K);
ber = zeros(3, n_chan, length(snr));

% params for frequency selective channel models
Ts = 1e-3;
Fd=0;
tau = [0 Ts Ts/2]; 
pdb = [0 -5 -10]; 

% number of fft/ifft points
npt = 64; 

% 802.11a frame params
N = 64; % number of subcarriers per ofdm burst
N_d_syms = 48; % number of data symbols per ofdm burst
N_gi = 12; % number of gaurd symbols per ofdm burst
mu = 16; % length of cyclic prefix
N_samp = N + mu; % total number of symbols per ofdm burst

% frame indicies 
pilot_idx = [11 25 39 53]; 
tmp = [6:58];
data_idx = setdiff(setdiff(tmp, pilot_idx), 27+5);

%% zf ofdm

for ii = 1:n_chan

    % generate and modulate random message
    bits = randi([0,M-1],1, N_syms*N_d_syms);
    mod_data = qammod(bits,M);  
    
    % set up ofdm frame
    ofdm_data = reshape(mod_data,N_d_syms,[]);
    ofdm_frame = zeros(N,N_syms);
    ofdm_frame(data_idx,:) = ofdm_data;
    ofdm_frame(pilot_idx,:) = 1;
    
    % add cyclic prefix 
    ifft_ofdm = ifft(ofdm_frame,npt);
    ofdm_tx = [ifft_ofdm(N_d_syms+1:N,:); ifft_ofdm];
      
    % generate random rayleigh channel
    Hchan = rayleighchan(Ts, Fd, tau, pdb); 
    chan = zeros(N_samp,N_syms);
    ofdm_tx_chan = zeros(N_samp,N_syms);
    for kk=1:N_syms
        chan(:,kk) = filter(Hchan,ones(N_samp,1));
        ofdm_tx_chan(:,kk) = chan(:,kk).*ofdm_tx(:,kk); %apply channel to signal
    end

    for jj = 1:length(snr)
        
        % add awgn noise
        n = sqrt(1/2)*(randn(N_samp, N_syms)...
            +1j*randn(N_samp, N_syms));
        
        rxchan = ofdm_tx_chan + 10^(-1*snr(jj)/20)*n;
        
        % deconstruct ofdm frame
        rx_rmgi = rxchan(mu+1:end,:);
        rx_ofdm_frame = fft(rx_rmgi, npt);
        
        % zf equilizer 
        rx_ofdm_frame_zf = rx_ofdm_frame ./ chan(mu+1:end,:);
        
        % un-pack data, demodulate
        rx_ofdm_data = rx_ofdm_frame_zf(data_idx,:);
        rx_ofdm_mod = reshape(rx_ofdm_data, 1, []);
        rx = qamdemod(rx_ofdm_mod, M);

        % compute ber curve
        [~, ber(1,ii,jj)] = biterr(bits, rx);

   end
end

%% mmse ofdm

for ii = 1:n_chan

    % generate and modulate random message
    bits = randi([0,M-1],1, N_syms*N_d_syms);
    mod_data = qammod(bits,M);  
    
    % set up ofdm frame
    ofdm_data = reshape(mod_data,N_d_syms,[]);
    ofdm_frame = zeros(N,N_syms);
    ofdm_frame(data_idx,:) = ofdm_data;
    ofdm_frame(pilot_idx,:) = 1;
    
    % add cyclic prefix 
    ifft_ofdm = ifft(ofdm_frame,npt);
    ofdm_tx = [ifft_ofdm(N_d_syms+1:N,:); ifft_ofdm];
      
    % select channel model
    Hchan = rayleighchan(Ts, Fd, tau, pdb);
    chan = zeros(80,N_syms);
    ofdm_tx_chan = zeros(80,N_syms);
    for kk=1:N_syms
        chan(:,kk) = filter(Hchan,ones(80,1));
        ofdm_tx_chan(:,kk) = chan(:,kk).*ofdm_tx(:,kk); %apply channel to signal
    end

    for jj = 1:length(snr)
        
        % add awgn noise
        n = sqrt(1/2)*(randn(N_samp, N_syms)...
            +1j*randn(N_samp, N_syms));
        
        rxchan = ofdm_tx_chan + 10^(-1*snr(jj)/20)*n;
        
        % deconstruct ofdm frame
        rx_rmgi = rxchan(mu+1:end,:);
        rx_ofdm_frame = fft(rx_rmgi, npt);
        
        % mmse equilizer 
        norm = conj(chan(mu+1:end,:)).*chan(mu+1:end,:) + ...
            10^(-1*snr(jj)/20);
        
        rx_ofdm_frame_zf = rx_ofdm_frame.*conj(chan(mu+1:end,:)) ./ norm;
        
        % un-pack data, demodulate
        rx_ofdm_data = rx_ofdm_frame_zf(data_idx,:);
        rx_ofdm_mod = reshape(rx_ofdm_data, 1, []);
        rx = qamdemod(rx_ofdm_mod, M);

        % compute ber curve
        [~, ber(2,ii,jj)] = biterr(bits, rx);

   end
end

%% plot

ber = mean(ber,2); %take mean across all channels

%plotting ber for different equalization techniques
semilogy(EbNo,ber(1,:),'DisplayName','zf')
hold on;
semilogy(EbNo,ber(2,:),'DisplayName','mmse');
hold on;

title('BER plots for OFDM equalization schemes');
xlabel('EbNo(dB)');
ylabel('Bit Error Rate');
legend('show');
