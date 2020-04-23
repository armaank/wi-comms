%% MIMO Link
clc; clear all;
%% simulation params

% modulation order
M = 16;
K = log2(M);     

% symbols per packet 
n_syms = 1e4;

% number of iterations per simulaton
n_iter = 3;

% 2x2 mimo link
Ntx = 2;
Nrx = 2;

EbNo = -5:2:25;                  
snr = EbNo + 10*log10(K);
ber = zeros(3, n_iter, length(snr));

% gaussian channel
H1 = sqrt(1/2)*(randn(Nrx, Ntx, n_syms*K, n_iter) +...
    1j*randn(Nrx, Ntx, n_syms*K, n_iter));

%% precoding mimo

U = zeros(Nrx, Ntx, n_syms*K);
S = zeros(Nrx, Ntx, n_syms*K);
V = zeros(Nrx, Nrx, n_syms*K);

precode = zeros(Nrx, 1, n_syms*K);
postcode = zeros(Nrx, 1, n_syms*K);

for ii = 1:n_iter

    % generate and modulate random message
    bits = randi([0 M-1], Ntx, 1, n_syms*K);     
    tx = qammod(bits, M); 
    
    % apply channel and perform pre-coding
    for kk=1:n_syms*K
        [U(:,:,kk),S(:,:,kk),V(:,:,kk)] = svd(H1(:,:,kk,ii));
        precode(:,:,kk) = H1(:,:,kk,ii)*V(:,:,kk)*tx(:,:,kk);
    end

    for jj = 1:length(snr)
        
        n = sqrt(1/2)*(randn(Nrx, 1, n_syms*K)...
            +1j*randn(Nrx, 1, n_syms*K));
        
        tx = precode + 10^(-1*snr(jj)/20)*n;
        
        % post-coding
        for kk=1:n_syms*K
            postcode(:,:,kk) = (S(:,:,kk)^-1)*U(:,:,kk)'*tx(:,:,kk);
        end
        
        rx = qamdemod(postcode,M);

        % compute ber curve
        [~, ber(1,ii,jj)] = biterr(bits, rx);

   end
end

%% zf mimo

zf_eq = zeros(Nrx, 1, n_syms*K);
txchan = zeros(Nrx, 1, n_syms*K); 
W = zeros(Nrx, Ntx, n_syms*K);

for ii = 1:n_iter

    % generate and modulate random message
    bits = randi([0 M-1], Ntx, 1, n_syms*K);     
    tx = qammod(bits, M); 
    
    % apply channel 
    for kk=1:n_syms*K
        txchan(:,:,kk) = H1(:,:,kk,ii)*tx(:,:,kk);
    end

    for jj = 1:length(snr)
        
        n = sqrt(1/2)*(randn(Nrx, 1, n_syms*K)...
            +1j*randn(Nrx, 1, n_syms*K));
        
        tx = txchan + 10^(-1*snr(jj)/20)*n;
        
        % zf equalizer
        for kk=1:n_syms*K
            W(:,:,kk) = (H1(:,:,kk,ii)'*H1(:,:,kk,ii))^-1*H1(:,:,kk,ii)';
            zf_eq(:,:,kk) = W(:,:,kk)*tx(:,:,kk);
        end
        rx = qamdemod(zf_eq,M);

        % compute ber curve
        [~, ber(2,ii,jj)] = biterr(bits, rx);

   end
end

%% mmse mimo


%% ploting

ber = mean(ber,2); %take mean across all iterations

%plotting ber for different equalization techniques
semilogy(EbNo,ber(1,:),'DisplayName','precoding');
hold on;
semilogy(EbNo,ber(2,:),'DisplayName','zero-forcing');
hold on;
semilogy(EbNo,ber(3,:),'DisplayName','mmse');

title('BER plots for MIMO schemes');
xlabel('EbNo (dB)');
ylabel('Bit Error Rate');
legend('show');
                  




