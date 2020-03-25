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
seed = ones(8,1);
M_seq = lfsr(seq_len, poly, seed);

% cfs for rrc filter
b_rrc = [0.0038; 0.0052; -0.0044; -0.0121; -0.0023; 0.0143; 0.0044;...
    -0.0385; -0.0563; 0.0363; 0.2554; 0.4968; 0.6025; 0.4968; .2554; ...
    0.0363; -0.0563; -0.0385; 0.0044; 0.0143; -0.0023; -0.0121; ...
    -0.0044; 0.0052; 0.0038]; 

% generating Walsh channel
H = hadamard(8);

%% cdma decoding 

% filter recieved signal 
rx_filter = filter(b_rrc,1, rx);
n_frames = length(rx)/(o_samp*cpf);

% decimate filtered signal by oversampling factor
rx_filter_dec = downsample(rx_filter,o_samp);

% using cross correlation of rx w/ m-sequence to determine offset
Exy = xcorr(abs(rx_filter_dec(1:cpf)), M_seq);
offset = find(Exy > .1 , 1, 'first');

% extracting, demodulating first and last pilots
pilot_first = rx_filter_dec(offset: cpf + offset -1);
pilot_last = rx_filter_dec(offset + (n_frames-1)*cpf:end);

pilot_first_demod = pskdemod(pilot_first.');
pilot_last_demod = pskdemod(pilot_last.');

% extracting phase to determine phase, frequency offsets
pilot_first_phase = unwrap(angle(pilot_first));
pilot_last_phase = unwrap(angle(pilot_last));

pilot_first_freq = pilot_first_phase.' - pilot_first_demod*pi;
pilot_last_freq = pilot_last_phase.' - pilot_last_demod*pi;

pilot_first_freq_diff = mod(diff(pilot_first_freq), 2*pi);
pilot_last_freq_diff = mod(diff(pilot_last_freq), 2*pi);
freq_offset_pc = median([pilot_first_freq_diff; pilot_last_freq_diff]);
freq_offset = cr* freq_offset_pc / (2*pi)
