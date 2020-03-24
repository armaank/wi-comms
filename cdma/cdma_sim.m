%% CDMA Receiver

clear all; clc;
%% simulation params

% load in received signal
rx = load('./data/Rcvd_Kohli.mat');
rx = rx.Rcvd; % received signal

% params
cr = 1e6; % code rate
o_samp = 4; % oversample
rrc_rolloff = .75; % roll off for RRC filter
M = 2; % bpsk modulaton

% instantiating modulators
pskmod = comm.PSKModulator(M,0);    
pskdemod = comm.PSKDemodulator(M,0);

% M sequence
seq_len = 255; 
poly = [8 7 6 1];
seed = ones(8,1);
M_seq = lfsr(seq_len, poly, seed);

% cfs for rrc filter
b_rrc = [0.0038; 0.0052; -0.0044; -0.0121; -0.0023; 0.0143; 0.0044;...
    -0.0385; -0.0563; 0.0363; 0.2554; 0.4968; 0.6025; 0.4968; .2554; ...
    0.0363; -0.0563; -0.0385; 0.0044; 0.0143; -0.0023; -0.0121; ...
    -0.0044; 0.0052; 0.0038]; 

% 
H = hadamard(8);

%% simulation

% filter recieved signal 
rx_filter = filter(b_rrc,1, rx);

% decimate filtered signal by oversampling factor
rx_filter_dec = downsample(rx_filter,o_samp);





