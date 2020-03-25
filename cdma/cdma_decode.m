%% CDMA Receiver

clear all; clc;
%% params

% load in received signal
rx = load('./data/Rcvd_Kohli.mat');
rx = rx.Rcvd; % received signal

% params
cr = 1e6; % chip rate
o_samp = 4; % oversample
rrc_rolloff = .75; % roll off for RRC filter
M = 2; % bpsk modulaton
cpf = 255; % chips per frame

% instantiating modulators
pskmod = comm.PSKModulator(M,0);    
pskdemod = comm.PSKDemodulator(M,0);

% M-sequence
seq_len = cpf; 
poly = [8 7 6 1];
seed = [1, zeros(1,7)];
M_seq = lfsr(seq_len, poly, seed);

% cfs for rrc filter
b_rrc = [0.0038; 0.0052; -0.0044; -0.0121; -0.0023; 0.0143; 0.0044;...
    -0.0385; -0.0563; 0.0363; 0.2554; 0.4968; 0.6025; 0.4968; .2554; ...
    0.0363; -0.0563; -0.0385; 0.0044; 0.0143; -0.0023; -0.0121; ...
    -0.0044; 0.0052; 0.0038]; 

% generating Walsh channel
H = hadamard(8);

% generating test pilot sequence
pilot = zeros(8*4,1);
pilot_mod = pskmod(pilot);
pilot_coded = (pilot_mod*H(1,:));
pilot_rx = xor(pskdemod(pilot_coded(1:end-1).'),M_seq).';

pskmod.release();
pskdemod.release()

pilot_rx_mod = -(pskmod(pilot_rx.')).';

%% cdma decoding 

% filter recieved signal 
rx_filter = filter(b_rrc,1, rx);
n_frames = length(rx)/(o_samp*cpf);

% decimate filtered signal by oversampling factor
rx_filter_dec = downsample(rx_filter,o_samp);

% using cross correlation of rx w/ m-sequence to determine offset
Exy = xcorr(abs(rx_filter_dec(1:cpf)), M_seq);
offset = find(Exy > 1 , 1, 'first');

% extracting, demodulating first and last pilots
pilot_first = rx_filter_dec(offset: cpf + offset -1);
pilot_last = rx_filter_dec(offset + (n_frames-1)*cpf:end);

pilot_first_demod = pskdemod(-pilot_first.').';
pilot_last_demod = pskdemod(-pilot_last.').';

% extracting phase to determine phase, frequency offsets
pilot_first_phase = unwrap(angle(pilot_first));
pilot_last_phase = unwrap(angle(pilot_last));

pilot_first_freq = pilot_first_phase - pilot_first_demod*pi;
pilot_last_freq = pilot_last_phase - pilot_last_demod*pi;

pilot_first_freq_diff = mod(diff(pilot_first_freq), 2*pi);
pilot_last_freq_diff = mod(diff(pilot_last_freq), 2*pi);
freq_offset_pc = median([pilot_first_freq_diff, pilot_last_freq_diff]);

freq_offset_Hz = cr* freq_offset_pc / (2*pi)

% re-shift the pilot in ordinance w/ the frequency offset to compute phase
pilot_first_freq_offset = pilot_first.*...
    exp(-1j*freq_offset_pc*(0:length(pilot_first)-1));

phase_offset = median(angle(pilot_first.*...
    exp(-1i*freq_offset_pc*(0:length(pilot_first)-1))));
 
phase_offset_Deg = phase_offset*180 / pi

% extract data, reverse frequency and phase offset
data_rx_offset = rx_filter_dec(cpf+offset : offset - 1 + (n_frames-1)*cpf);

data_rx = data_rx_offset...
    .*exp(-1j*freq_offset_pc*(cpf:cpf+length(data_rx_offset)-1))...
     .*exp(-1j*phase_offset);

% despread and decode data
data_pn = data_rx - repmat(pilot_rx_mod,1,n_frames-2);
data_pn = bpskcdma(data_pn);

frames = reshape(data_pn, cpf, []);

data_frames_pn = frames(:,1:end-1);
carrier_frame_pn = frames(:,end); 
data_frames = zeros(192,size(data_frames_pn,2));

% demodulate and unscramble data and carrier
for ii = 1:(8*4)
    
    data_frame_demod = data_frames_pn(find(data_frames_pn(:,ii)),ii)>0;
    
    data_frames(:,ii) = xor(data_frame_demod,M_seq(1:192));

end

carrier_frame_demod = pskdemod(-carrier_frame_pn(find(carrier_frame_pn)));
carrier_frame = xor(carrier_frame_demod,...
    M_seq(1:length(carrier_frame_demod)));

% decode final frames
data = [data_frames(:);carrier_frame(:)];

pskmod.release()
pskdemod.release()

bits = pskdemod((-H(6,:)*reshape(-pskmod((data)).',8,[])/8).');
msg = char(bi2de(reshape(bits,8,[]).').')