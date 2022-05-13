% Universal constants used
k_b = 1.380649*10e-23; % m^2 kg s^-1 K^-1, Boltzmann constant
c0 = 3e8;        % m/s, speed of light
pi = 3.1415;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Antenna User input
SNRdB_user_input = (190);
SNR_dB = zeros(1,length(SNRdB_user_input)); % dB
counter= 1;
for v = SNRdB_user_input  
    SNR_dB(1,counter) = v;
    counter= counter+ 1;
end
P_tx_sat = 80 ; %Watts, power of the signal transmitted by satellite
P_tx_sat_dBW = 10*log10(P_tx_sat); % dBW, transmitted power in dBW
P_tx_sat_dBm = P_tx_sat_dBW + 30; % dBm, transmitted power in dBm

g_tx_dB = 17.8; % dBi , Gain of transmit antenna
g_rx_dB = 20.0; % dBi , Gain of receive antenna

fc = 20e9;      % Hz, carrier frequency 20GHz
B_c = 2e9; % Hz, carrier Bandwidth (24-26GHz)
lambda = c0/fc;  % m,  wavelength
nu  = 2*pi*fc/c0 ; % wavenumber of carrier signal
Temp = 70; % Kelvin, system temperature

P_noise_dbW = -120 ; % dBW , noise power
P_noise_dBm = P_noise_dbW + 30; % dBm, noise power
B_n = 30e6 ; % hz , Noise Bandwidth
sigma_sq = k_b * Temp * B_n ; % variance of AWGN


% 01/17/2022 10:20:36 UTC (MM/DD/YYYY HH:MM:SS)
utc = [2022 1 17 10 20 36];

% Receiver Antenna locations in Lat., Lon., Alt.
rcvr_1=[53.10657108482692, 8.850392209543823 18]; % Bremen Uni 53°N , 8.8°E, 18m 
Nrx= 8; % Create Nrx X Nrx MIMO at Ground station
gap = 1; % meters, gap between the Rx antennas one half wavelength


% satellite positions
inter_satellite_distance = (50:100:100000); % meters





