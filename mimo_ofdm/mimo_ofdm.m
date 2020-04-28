%% MIMO Link w/ OFDM
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
snr = EbNo + 10*log10(K) + 10*log10(64/80);
ber = zeros(2, 3, n_chan, length(snr));

% params for frequency selective channel models
Ts = 1e-3;
Fd = 0;
tau = [0 Ts Ts/2]; 
pdb = [0 -5 -10]; 

% 2x2 mimo link
Ntx = 2;
Nrx = 2;

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
data_idx = [1:5 7:19 21:26 28:33 35:47 49:53]+5; %frame parameters

H = sqrt(1/2)*(randn(Nrx, Ntx, N_syms*N_d_syms/2, n_chan) +...
    1j*randn(Nrx, Ntx, N_syms*N_d_syms/2, n_chan));

%% ofdm zero forcing

U = zeros(Nrx, Ntx, N_syms*N_d_syms/2);
S = zeros(Nrx, Ntx, N_syms*N_d_syms/2);
V = zeros(Nrx, Ntx, N_syms*N_d_syms/2);

precode = zeros(Nrx,1,N_syms*N_d_syms/2);
txchan = zeros(Nrx,1,N_syms*N_d_syms/2);
postcode = zeros(Nrx,1,N_syms*N_d_syms/2);

W = zeros(Nrx, Ntx, N_syms*N_d_syms/2);
zf_eq = zeros(Nrx, 1, N_syms*N_d_syms/2);
mmse_eq = zeros(Nrx, 1, N_syms*N_d_syms/2);


for ii = 1:n_chan
    
    % generate and modulate random message
    bits = randi([0,M-1],1, N_syms*N_d_syms);     
    bits_comp = reshape(bits,2,[]);
    mod_data = qammod(bits,M);
    
    % set up ofdm frame
    ofdm_data = reshape(mod_data,N_d_syms,[]);
    ofdm_frame = zeros(N,N_syms);
    ofdm_frame(data_idx,:) = ofdm_data; 
    ofdm_frame(pilot_idx,:) = 1;
    
    % add cyclic prefix
    ifft_ofdm = ifft(ofdm_frame,64);
    ofdm_trans = [ifft_ofdm(N_d_syms+1:N,:); ifft_ofdm];

    % generate random rayleigh channel
    Hchan = rayleighchan(Ts, Fd, tau, pdb);
    chan = zeros(N_samp,N_syms);
    ofdm_tx_chan = zeros(N_samp,N_syms);
    for kk=1:N_syms
        chan(:,kk) = filter(Hchan,ones(N_samp,1));
        ofdm_tx_chan(:,kk) = chan(:,kk).*ofdm_trans(:,kk); 
    end
    
    for jj = 1:length(snr)
        
        % add awgn noise
        n = sqrt(1/2)*(randn(N_samp,N_syms)+1j*randn(N_samp,N_syms));
        
        rxchan = ofdm_tx_chan + 10^(-1*(snr(jj)-3)/20)*n;

        % deconstruct ofdm frame
        rx_rmgi = rxchan(mu+1:end,:);
        rx_ofdm_frame = fft(rx_rmgi,N);
        
        % zf eq for ofdm
        rx_ofdm_frame_zf = rx_ofdm_frame./chan(mu+1:end,:);
        
        % unpack data
        rx_ofdm_data = rx_ofdm_frame_zf(data_idx,:);
        rx_ofdm_mod = reshape(rx_ofdm_data,1,[]); 
        
        % mimo link
        tx_mimo = reshape(rx_ofdm_mod,2,[]);
        tx_mimo = permute(tx_mimo,[1,3,2]);
        
        % apply channel and perform precoding 
        for kk=1:N_syms*N_d_syms/2 
            
            [U(:,:,kk),S(:,:,kk),V(:,:,kk)] = svd(H(:,:,kk,ii));
            precode(:,:,kk) = H(:,:,kk,ii)*V(:,:,kk)*tx_mimo(:,:,kk);
            
        end
        
        n = sqrt(1/2)*(randn(Nrx,1,N_syms*N_d_syms/2)...
            +1j*randn(Nrx,1,N_syms*N_d_syms/2));
        
        tx = precode + 10^(-1*(snr(jj)-3)/20)*n;
        
        % post-coding
        for kk=1:N_syms*N_d_syms/2
            
            postcode(:,:,kk) = S(:,:,kk)^-1*U(:,:,kk)'*tx(:,:,kk);
            
        end

        rx = qamdemod(postcode,M);
        % compute ber curve for mimo precoding 
        [~, ber(1,1,ii,jj)] = biterr(bits_comp, squeeze(rx)); 
        
        % apply channel (mimo zf) 
        for kk=1:N_syms*N_d_syms/2
            txchan(:,:,kk) = H(:,:,kk,ii)*tx_mimo(:,:,kk);
        end
        
        tx = txchan + 10^(-1*(snr(jj)-3)/20)*n;

        % zf equalizer 
        for kk=1:N_syms*N_d_syms/2
            
            W(:,:,kk) = (H(:,:,kk,ii)'*H(:,:,kk,ii))^-1*H(:,:,kk,ii)';
            zf_eq(:,:,kk) = W(:,:,kk)*tx(:,:,kk);
            
        end
        
        rx = qamdemod(zf_eq,M);
        
        % compute ber curve for mimo zf eq.
        [~, ber(1,2,ii,jj)] = biterr(bits_comp, squeeze(rx)); 

        % mmse equalizer 
        for kk=1:N_syms*N_d_syms/2
            
            W(:,:,kk) = (H(:,:,kk,ii)'*H(:,:,kk,ii)+...
                eye(Ntx)*10^(-1*(snr(jj)-3)/20))^-1*H(:,:,kk,ii)';
            
            mmse_eq(:,:,kk) = W(:,:,kk)*tx(:,:,kk);
            
        end
        
        rx = qamdemod(mmse_eq,M);

        % Compute and store the BER for mimo mmse eq. 
        [~, ber(1,3,ii,jj)] = biterr(bits_comp, squeeze(rx)); 
        
    end
end



%% mmse OFDM

%% plot

ber = mean(ber, n_chan); %take mean across all channels

%plotting 
figure;
semilogy(EbNo,squeeze(ber(1,1,1,:)),'DisplayName','precoding');
hold on;
semilogy(EbNo,squeeze(ber(1,2,1,:)),'DisplayName','zf');
semilogy(EbNo,squeeze(ber(1,3,1,:)),'DisplayName','mmse');

title('BER plots for MIMO Schemes w/ OFDM zf');
xlabel('EbNo(dB)');
ylabel('Bit Error Rate');
legend('show');
hold off;

