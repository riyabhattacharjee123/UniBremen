clear;
close all;
clear;
clc;
%run('User_input_data.m');
%run('mimo3_16.m');
%run('mimoX_16.m');


run('User_input_data.m');
BER_plot = zeros(1,length(SNRdB_user_input));
R_sum_rate_plot = zeros(1,length(SNRdB_user_input));
R_sum_rate_plot = zeros(1,length(theta_el_deg_user));


for te = 1:length(theta_el_deg_user)
    
    theta_el_deg = theta_el_deg_user(te)     
    run('get_eci_loc_sat.m');
    run('receiver_antenna_locations.m');

    %inter_sat_dist = 2000 ; % meters distance_1 50:10:100000 100:1:100

    %side = inter_sat_dist
    %run('GSfixedcoordinates.m');
    %run('cartesianTOspherical.m'); 
        %run('steeringVectors.m');
        %run('channel_rate3Dmax.m');       
    %run('calc_inter_satellite_distance_normal.m');


    % satellite positions
    %run('calc_inter_sat_distance.m'); % based on our condition
    %inter_satellite_distance = 50:300:1000000 ; % meters distance_1 50:10:100000 100:1:100
    
    for isd = 1: length(inter_satellite_distance)
        side = inter_satellite_distance(isd)
        for su = 1:length(SNRdB_user_input)

            SNR_dB(1,su) = SNRdB_user_input(1,su);    
            SNR_lin(1,su) = 10^(SNRdB_user_input(1,su)/10);
            SNR_linear = SNR_lin(1,su);

            PSD_N_0 = P_tx_gs / SNR_linear ; % Power Spectral Density of noise

            run('GSfixedcoordinates.m');
            run('cartesianTOspherical_r_start_sat.m'); 
            run('steeringVectors.m');
            %run('channel_rate3Dmax.m');  

            %Y1(isd)= real(R_opt) ; %gs_sat_polar_theta(1,1) * (180/pi) ;
            X1(isd)= inter_satellite_distance(isd)/1000;  
            %Y2(isd)= real(sum_dft_1) ;    
            run('linear_receiver.m'); 
            run('linear_receiver_sinr.m');
            run('linear_receiver_ber.m');
            BER_plot(su)= BER(1,1);
            R_sum_rate_plot(te)= R_sum_rate;             

            run('LMS_adaptive_filter.m'); 
        end

    end
end

% plot graph
figure();
plot(X1,real(R_sum_rate_plot),'*-r',X1,imag(R_sum_rate_plot),'*-g');
hold on;
ylabel("Sum Rate (bps/Hz)");
xlabel("Inter satellite distance (Km)");
title("Distance between satellites vs Sum Rate");
legend('Re(3 sat. Triangle form.)','Im(3 sat. Triangle form.)')
grid on;


%figure();
%plot(X1,Y2);
%ylabel("DFT matrix sum");
%xlabel("Inter satellite distance (Km)");
%title("Distance between satellites vs sum of DFT matrix");
%legend('3 X 64 MIMO')
%grid on;

figure();
semilogy(SNR_dB,BER_plot,'*-');
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('SNR Vs BER plot for QPSK Modualtion in Rayleigh Channel for satellite 1 data (3 sat. Triangular formation)');

figure();
plot(SNR_dB,real(R_sum_rate_plot),'*-r',SNR_dB,imag(R_sum_rate_plot),'*-g');
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Sum Rate (bps/Hz)');
legend('Re(Sum Rate)','Im(Sum Rate)');
title('SNR Vs Sum Rate plot for QPSK Modualtion in Rayleigh Channel for satellite 1 data (3 sat. Triangular formation)');

figure();
plot(theta_el_deg_user,real(R_sum_rate_plot),'*-r',theta_el_deg_user,imag(R_sum_rate_plot),'*-b');
hold on;
ylabel("Sum Rate (bps/Hz)");
xlabel("Elevation angle (deg)");
title("Elevation angle vs Sum Rate");
legend('Re(3 sat. Triangle form.)','Im(3 sat. Triangle form.)')
grid on;

