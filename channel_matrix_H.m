clear all
close all
clc
% channel_matrix_H.m
% we calculate channel matrix [H] (NtxNr)
% Number of transmitting antennas (Nr) = 3
% Number of receivng antennas (Nr) = 3
% Hence [H] is a 3X3 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal constants used
k_b = 1.380649*10e-23; % m^2 kg s^-1 K^-1, Boltzmann constant
c0 = 3e8;        % m/s, speed of light
pi = 22/7;

%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\coordinates_convertor.m');
run('coordinates_convertor.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Antenna User input variables

Nr = 3; % number of receiver antenna
Nt = 3; % number of transmitting antenna

SNR_dB = zeros(1,41); % dB
counter= 1;
for v = 50:10:450
   SNR_dB(1,counter) = v;
   counter= counter+ 1;
end

P_tx_sat = 180 ; %Watts, power of the signal transmitted by satellite
P_tx_sat_dBW = 10*log10(P_tx_sat); % dBW, transmitted power in dBW
P_tx_sat_dBm = P_tx_sat_dBW + 30; % dBm, transmitted power in dBm

g_tx_dB = 17.8; % dBi , Gain of transmit antenna
g_rx_dB = 20.0; % dBi , Gain of receive antenna

fc = 25e9;      % Hz, carrier frequency 25GHz
B_c = 2e9; % Hz, carrier Bandwidth (24-26GHz)
lambda = c0/fc;  % m,  wavelength
nu  = 2*pi*fc/c0 ; % wavenumber of carrier signal
Temp = 70; % Kelvin, system temperature

P_noise_dbW = -120 ; % dBW , noise power
P_noise_dBm = P_noise_dbW + 30; % dBm, noise power
B_n = 30e6 ; % hz , Noise Bandwidth
sigma_sq = k_b * Temp * B_n ; % variance of AWGN
%snr_rx = P_tx_sat/(10^(P_noise_dbW/10));
%snr_rx_db = 10*log10(snr_rx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Channel matrix and SNR matrix for each channel
H = [];
for n=1:size(rcvr_pos_eci,1)
    for m=1:size(trx_pos_eci,1)
               
        del_x = trx_pos_eci(m,1)- rcvr_pos_eci(n,1); % meters
        del_y = trx_pos_eci(m,2)- rcvr_pos_eci(n,2); % meters
        del_z = trx_pos_eci(m,3)- rcvr_pos_eci(n,3); % meters        
        distance(n,m) = sqrt(del_x^2+del_y^2+del_z^2); % meters
        
        % calculate the free space loss (fsl)
        fsl_dB(n,m) = -147.55 + 20*log10(distance(n,m)) + 20*log10(fc); % dB
        fsl_linear(n,m) = 10^(fsl_dB(n,m)/10); % linear scale
        
        % calculate the pathloss in dB
        pathloss_dB(n,m) = -20*log10(lambda/(4*pi*distance(n,m)))- g_tx_dB - g_rx_dB;        
        pathloss_linear(n,m) = 10^(pathloss_dB(n,m)/10);
        
        % calculate the elements of channel matrix
        H(n,m) = (1/sqrt(pathloss_linear(n,m)))* exp(-1i*nu*distance(n,m));
       % H(n,m) = (1/sqrt(fsl_linear(n,m)))* exp(-1i*nu*distance);
        
        % calculate elements of SNR for each link
        SNR_channel(n,m) = P_tx_sat*pathloss_linear(n,m)/sigma_sq ;
        SNR_channel_dB(n,m)= 10*log10(SNR_channel(n,m));
        
    end    
end

% Signal transmitted
M = 8; % number of bits in a message stream
k = log2(M);% Number of bits per symbol
x_signal = randi([0,1],Nt,M); % generating random binary stream of M bits

% generate AWGN noise
% convert Noise power from dBW to linear scale
P_noise_linear = 10^(P_noise_dbW / 10); % Watts , sigma^2
for p=1:Nr
    for q=1:M
        noise_signal(p,q) = P_noise_linear;
    end
end

% received signal
y_signal = H*x_signal + noise_signal; %Received Signal at ground station
disp(y_signal);

% Channel Rate
H_hermitian = transpose(conj(H));
I = ones(Nr,Nr); % Identity matrix
% Calculate the Channel capacity for each SNR
for s = 1:size(SNR_dB,2)  
    SNR = 10^(SNR_dB(s)/10); % converting SNR from dB to linear scale
    R =  log2(det(I+((SNR/Nt)*H*H_hermitian))); % bits/sec/Hz
    
    X1(s)=SNR;
    X2(s)= SNR_dB(s);
    Y(s)=R;
end   

% plot graph
figure()
plot(X1,Y);
xlabel("SNR in linear scale");
ylabel("Channel Rate (R) (bps/Hz)");
title("Channel rate vs SNR");
legend('3X3 MIMO')
grid on;

figure();
plot(X2,Y);
xlabel("SNR in dB");
ylabel("Channel Rate (R) (bps/Hz)");
title("Channel rate vs SNR_dB");
legend('3X3 MIMO')
grid on;





