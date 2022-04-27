clear all
% channel_matrix_H.m
% we calculate channel matrix [H] (NtxNr)
% Number of transmitting antennas (Nr) = 3
% Number of receivng antennas (Nr) = 3
% Hence [H] is a 3X3 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal constants used
k_b = 1.380649*10e-23; % m^2 kg s^-1 K^-1, Boltzmann constant
c0 = 3e8;        % m/s, propagation speed

run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\coordinates_convertor.m');

% Antenna parameters user input
Nr = 3; % number of receiver antenna
Nt = 3; % number of transmitting antenna

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Antenna User input variables

SNR_dB = zeros(1,20); % dB
counter= 1;
for v = 10:5:100
   SNR_dB(1,counter) = v;
   counter= counter+ 1;
end

P_noise_dbW = -120 ; % dBW , noise power
%P_tx_dBW = SNR_dBW + P_noise_dbW ; % dBW, transmit signal power
P_tx_sat = 250 ; %Watts, power of the signal transmitted by satellite
g_tx_dB = 17.8; % dBi , Gain of transmit antenna
g_rx_dB = 20.0; % dBi , Gain of receive antenna

fc = 25e9;      % Hz, carrier frequency 24GHz
B_c = 2e9; % Hz, carrier Bandwidth (24-26GHz)
lambda = c0/fc;  % m,  wavelength
nu  = 2*pi*fc/c0 ; % wavenumber of carrier signal
Temp = 70; % Kelvin, system temperature
B_n = 30e6 ; % hz , Noise Bandwidth
sigma_sq = k_b * Temp * B_n ; % variance of AWGN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Channel matrix and SNR matrix
H = [];
for n=1:size(rcvr_pos_eci,1)
    for m=1:size(trx_pos_eci,1)
               
        del_x = trx_pos_eci(m,1)- rcvr_pos_eci(n,1); % meters
        del_y = trx_pos_eci(m,2)- rcvr_pos_eci(n,2); % meters
        del_z = trx_pos_eci(m,3)- rcvr_pos_eci(n,3); % meters        
        distance = sqrt(del_x^2+del_y^2+del_z^2); % meters
        
        % calculate the pathloss in dB
        pathloss_dB = g_tx_dB+ g_rx_dB -20*log10(lambda/(4*pi*distance)); 
        
        pathloss_linear = 10^(pathloss_dB/10);
        
        % calculate the elements of channel matrix
        H(n,m) = (1/sqrt(pathloss_linear))* exp(-j*nu*distance);
        
        % calculate elements of SNR for each link
        %SNR(n,m) = P_tx_sat*pathloss_linear/sigma_sq ;
        
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
I = ones(Nr,Nr);
% Calculate the Channel capacity for each SNR
for s = 1:size(SNR_dB,2)    
    SNR = 10^(SNR_dB(s)/10); % converting SNR from dB to linear scale
    R =  log2(det(I+(SNR*H*H_hermitian))); % bits/sec/Hz
    X(s)=SNR
    Y(s)=R      
end   

% plot graph
plot(X,Y);
xlabel("SNR");
ylabel("Channel Rate (R) (bps/Hz)");
title("Channel rate vs SNR");
legend('3X3 MIMO')
grid on;






