% Universal constants used
k_b = 1.380649*10e-23; % m^2 kg s^-1 K^-1, Boltzmann constant
c0 = 3e8;        % m/s, speed of light
pi = 3.1415;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Antenna User input
SNRdB_user_input = (200); %dB
SNR_lin = 10^(SNRdB_user_input/10);
SNR_dB = zeros(1,length(SNRdB_user_input)); % dB
counter= 1;
for v = SNRdB_user_input  
    SNR_dB(1,counter) = v;
    counter= counter+ 1;
end
P_tx_sat = 100 ; %Watts, power of the signal transmitted by satellite
P_tx_sat_dBW = 10*log10(P_tx_sat); % dBW, transmitted power in dBW
P_tx_sat_dBm = P_tx_sat_dBW + 30; % dBm, transmitted power in dBm

g_tx_dB = 17.8; % dBi , Gain of transmit antenna
g_rx_dB = 20.0; % dBi , Gain of receive antenna

gain_channel_ul =  18.00; %-100.0 ; % dB , channel gain uplink
gain_channel_dl = -105.0 ; % dB , channel gain downlink

P_tx_gs_dBm = P_tx_sat_dBm ; % dBm 
P_tx_gs_dBW = P_tx_gs_dBm - 30; % dbW
P_tx_gs = 10^(P_tx_gs_dBW/10); % Watts

fc = 6e9; %20e9;      % Hz, carrier frequency 6GHz
B_c = 2e9; % Hz, carrier Bandwidth (24-26GHz)
lambda = c0/fc;  % m,  wavelength
nu  = 2*pi*fc/c0 ; % wavenumber of carrier signal
Temp = 70; % Kelvin, system temperature

P_noise_dbW = -200; % -135 ; % dBW , noise power
P_noise_dBm = P_noise_dbW + 30; % dBm, noise power
B_n = 30e6 ; % hz , Noise Bandwidth
P_noise_var = 10^(P_noise_dbW/100);
sigma_sq = k_b * Temp * B_n ; % variance of AWGN

% PSD_N_0 = P_tx_gs/SNR_lin
PSD_N_0 = P_tx_gs / SNR_lin ; % Power Spectral Density of noise



% 01/17/2022 10:20:36 UTC (MM/DD/YYYY HH:MM:SS)
utc = [2022 1 17 10 20 36];

% Receiver Antenna locations in Lat., Lon., Alt.
rcvr_1=[53.10657108482692, 8.850392209543823 18]; % Bremen Uni 53°N , 8.8°E, 18m 
%rcvr_1=[0.009813 -78.452791 18]; % Quito Ecuador 0°N , -78.8°W, 18m 
Nrx= 8; % Create Nrx X Nrx MIMO at Ground station
D_Ant = c0/(2*fc); % 1; % meters, gap between the Rx antennas one half wavelength, D_A


inter_satellite_distance = 400000 ; % meters 50:300:1000000

% Number of antenna attached to a satellite
N_sat_ant = 1 ;
% Signal transmitted
%M = 5; % number of bits in a message stream
%k = log2(M);% Number of bits per symbol
%x_signal = randi([0,1],Nrx*Nrx,M); % generating random binary stream of M bits

%t = 0.0001:0.0001:1;





