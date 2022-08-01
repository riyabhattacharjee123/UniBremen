clear;
close all;
clear;
clc;



run('User_input_data.m');
run('receiver_antenna_locations.m')

% 3 satellites
% satellite positions
%run('calc_inter_satellite_distance_normal.m'); % based on our condition
inter_satellite_distance = 50:1000:10000000; %50:10:50000  ; % meters distance_1 50:10:50000

for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors.m');
    run('channel_rate3Dmax.m');    
    
    Y1_3(isd)= real(R_opt) ; %gs_sat_polar_theta(1,1) * (180/pi) ;
    X1_3(isd)= inter_satellite_distance(isd)/1000;  
    %Y2(isd)= real(sum_dft_1) ;
    
end




% for 4 satellites
% satellite positions
run('calc_inter_sat_distance_4.m'); % based on our condition
inter_satellite_distance = 50:1000:10000000; %50:10:50000  ; % meters distance_1 50:10:100000

for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('GSfixedcoordinates_4.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors_4.m');
    run('channel_rate3Dmax.m');    
    
    Y1_4(isd)= real(R_opt) ; %gs_sat_polar_theta(1,1) * (180/pi) ;
    X1_4(isd)= inter_satellite_distance(isd)/1000;  
    %Y2(isd)= real(sum_dft_1) ;
    
end

% plot graph
figure();
%
plot(X1_3,Y1_3);
hold on;
plot(X1_4,Y1_4);
ylabel("Optimum Channel Rate (bps/Hz)");
xlabel("Inter satellite distance (Km)");
title("Distance between satellites vs R optimum");
legend('3 X 64 MIMO','4 X 64 MIMO')
%axis equal;
grid on;