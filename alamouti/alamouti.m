% space-time diversity via Alamouti Codes
% wireless link simulation over rayleigh flat-fading channel
clc; clear all;
%% simulation parameters

% modulation order
M = 2;

% bits per symbol
k = log2(M);

% instantiating modulators
pskmod = comm.PSKModulator(M,0);    
pskdemod = comm.PSKDemodulator(M,0);

EbNo = 0:2:50;                  
snr = EbNo + 10*log10(k);       

% messsage length 
N = 1e6; 

% generating Rayleigh flat-fading channels
fd = 10; 
raychan1 = rayleigh(fd, N);
raychan2 = rayleigh(fd, N);
raychan3 = rayleigh(fd, N);
raychan4 = rayleigh(fd, N);

%% simulation 1 - uncoded bpsk

% init 
v_demod = zeros(N,length(snr));
v_rx = zeros(N,length(snr));

% generate random message
v = randi([0 M-1],N,1);

% modulate
v_mod = step(pskmod,v);

% transmit through channel 
v_raychan = raychan1.*v_mod;

for ii=1:length(snr)
    
   v_rx(:,ii) = awgn(v_raychan,snr(ii),'measured')./raychan1;
   
   v_demod(:,ii) = step(pskdemod,v_rx(:,ii)); 
   
end
% compute ber
[~,ber1] = biterr(v_demod,v);

%% simulation 2 - bpsk w/ mrrc (2 rx)

% init
combiner = zeros(N,length(snr));
v_demod = zeros(N,length(snr));

% generate random message
v = randi([0 M-1],N,1);

% modulate
v_mod = step(pskmod,v);
% transmit through channel 

h = raychan1; 
v_raychan = [raychan1.*v_mod, raychan2.*v_mod]; 

for ii=1:length(snr)
    
   v_rx = awgn(v_raychan,snr(ii),'measured'); 
   
   % mrrc
   combiner(:,ii) = sum(conj(h).*v_rx,2)./sum(h.*conj(h),2);
   
   % demodulate
   v_demod(:,ii) = step(pskdemod,combiner(:,ii)); 
   
end
% compute ber
[~,ber2] = biterr(v_demod,v);

%% simulation 3 - bpsk w/ mrrc (4 rx)

% init
combiner = zeros(N,length(snr));
v_demod = zeros(N,length(snr));

% generate random message
v = randi([0 M-1],N,1);

% modulate
v_mod = step(pskmod,v);

% transmit through channel
v_raychan1 = raychan1.*v_mod;
v_raychan2 = raychan2.*v_mod; 
v_raychan3 = raychan3.*v_mod; 
v_raychan4 = raychan4.*v_mod; 

h = [raychan1,raychan2, raychan3, raychan4]; 
v_raychan = [v_raychan1, v_raychan2, v_raychan3, v_raychan4]; 

for ii=1:length(snr)
    
   v_rx = awgn(v_raychan,snr(ii),'measured'); 
   
   % mrrc
   combiner(:,ii) = sum(conj(h).*v_rx,2)./sum(h.*conj(h),2);
   
   % demodulate
   v_demod(:,ii) = step(pskdemod,combiner(:,ii)); 
   
end
% compute ber
[~,ber3] = biterr(v_demod,v);

%% simulation 4 - bpsk w/ Alamouti coding (2 tx, 1 rx)

% init
combiner = zeros(N,length(snr));
v_decoded = zeros(N,length(snr));
v_demod = zeros(N,length(snr));

% generate random message
v = randi([0 M-1],N,1);

% modulate
v_mod = step(pskmod,v);

% alamouti coding
v_coded=zeros(N,2);

s1=v_mod(1:2:end); 
s2=v_mod(2:2:end);

v_coded(1:2:end,:)=sqrt(0.5)*[s1,s2];
v_coded(2:2:end,:)=sqrt(0.5)*[-conj(s2),conj(s1)];

% transmit
h = sqrt(1/2)*kron([raychan1(1:N/2), raychan2(1:N/2)], [1;1]); 
v_raychan = h.*v_coded;

for ii=1:length(snr)

   v_rx = awgn(v_raychan,snr(ii),'measured'); 
   
   combiner(:,ii) = sum(v_rx,2);
   
   % alamouti decoding
   u1 = combiner(1:2:end,ii); u2 = combiner(2:2:end,ii);
   u = [kron(u1,[1;1]), kron(conj(u2),[1;1])]; 
   
   h_rx = zeros(N,2); 
   h1=h(1:2:end,1); 
   h2=h(1:2:end,2);
   h_rx(1:2:end,:) = [conj(h1), h2];
   h_rx(2:2:end,:) = [conj(h2), -h1];
   
   v_decoded(:,ii) = sum(h_rx.*u,2)./sum(h_rx.*conj(h_rx),2);
   
   % demodulate
   v_demod(:,ii) = step(pskdemod,v_decoded(:,ii)); 
   
end
% compute ber
[~,ber4] = biterr(v_demod,v);
semilogy(snr, ber4)
%% Simulate BPSK through rayleigh channel w/ Alamouti (2 Tx, 2 Rx)

% init
combiner = zeros(N,2,length(snr));
v_decoded = zeros(N,length(snr));
v_demod = zeros(N,length(snr));

% generate random message
v = randi([0 M-1],N,1);

% modulate
v_mod = step(pskmod,v);

% alamouti coding 
v_coded=zeros(N,4);

s1=v_mod(1:2:end); 
s2=v_mod(2:2:end);

v_coded(1:2:end,:)=sqrt(1/2)*[s1,s2, s1,s2];
v_coded(2:2:end,:)=sqrt(1/2)*[-conj(s2),conj(s1), -conj(s2),conj(s1)];

% transmit
h = sqrt(.5)*kron([raychan1(1:N/2), raychan2(1:N/2),...
    raychan3(1:N/2), raychan4(1:N/2)], [1;1]);
v_raychan = h.*v_coded;

for ii=1:length(snr)
    
   v_rx = awgn(v_raychan,snr(ii),'measured'); 
   
   % combiner
   combiner(:,:,ii) = [sum(v_rx(:,1:2),2), sum(v_rx(:,3:4),2)];
   
   % alamouti decoding
   u11 = combiner(1:2:end,1,ii); 
   u21 = combiner(2:2:end,1,ii);
   u12 = combiner(1:2:end,2,ii); 
   u22 = combiner(2:2:end,2,ii);
   
   u = [kron(u11,[1;1]),kron(conj(u21),[1;1]),...
       kron(u12,[1;1]),kron(conj(u22),[1;1])]; 
   
   h_rx = zeros(N,4); 
   h1=h(1:2:end,1); 
   h2=h(1:2:end,2); 
   h3=h(1:2:end,3); 
   h4=h(1:2:end,4);
   
   h_rx(1:2:end,:) = [conj(h1),h2,  conj(h3),h4];
   h_rx(2:2:end,:) = [conj(h2),-h1, conj(h4),-h3];
   
   v_decoded(:,ii) = sum(h_rx.*u,2)./sum(h_rx.*conj(h_rx),2);
   
   % demodulate
   v_demod(:,ii) = step(pskdemod,v_decoded(:,ii)); 
end
% compute ber
[~,ber5] = biterr(v_demod,v);

%% Plot

semilogy(snr,ber1,'-o',snr,ber2,'-v',snr,ber3,'-s',...
        snr,ber4,'-d',snr,ber5,'-^', 'LineWidth', 2);
title('BER for BPSK through  a Rayleigh Channel','FontSize', 14 );
grid on;
xlabel('Eb/No (dB)', 'FontSize', 14);
ylabel('Bit Error Rate', 'FontSize', 14);
legend('uncoded (1 Tx, 1 Rx)','MRRC (1 Tx, 2 Rx)',...
      'MRRC (1 Tx, 4 Rx)', 'Alamouti (2 Tx, 1 Rx)',...
      'Alamouti (2 Tx, 2 Rx)');
ax = gca;
ax.LineWidth = 1.75;






