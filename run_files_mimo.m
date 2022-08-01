clear;
close all;
clear;
clc;
%run('User_input_data.m');
%run('mimo3_16.m');
%run('mimoX_16.m');


run('User_input_data.m');
run('receiver_antenna_locations.m')

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
    run('GSfixedcoordinates.m');
    run('cartesianTOspherical_r_start_sat.m'); 
    run('steeringVectors.m');
    run('channel_rate3Dmax.m');  
    
    Y1(isd)= real(R_opt) ; %gs_sat_polar_theta(1,1) * (180/pi) ;
    X1(isd)= inter_satellite_distance(isd)/1000;  
    %Y2(isd)= real(sum_dft_1) ;
    
end

% plot graph
figure();
plot(X1,Y1);
ylabel("Optimum Channel Rate (bps/Hz)");
xlabel("Inter satellite distance (Km)");
title("Distance between satellites vs R optimum");
legend('3 sat. Triangle form.')
grid on;


%figure();
%plot(X1,Y2);
%ylabel("DFT matrix sum");
%xlabel("Inter satellite distance (Km)");
%title("Distance between satellites vs sum of DFT matrix");
%legend('3 X 64 MIMO')
%grid on;


