% Antenna locations
run('receiver_antenna_locations.m')

% satellite positions
inter_satellite_distance = (50:100:50000); % meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leader follower configuration
r_start_all_sat = []; % Initialize to null.
v_start_all_sat = []; % Initialize to null.
for ss = 1:length(inter_satellite_distance)
   side = inter_satellite_distance(ss)  
   run('gco_hcw_leader_follower.m');
   run('init_start.m');
   r_start_all_sat = [r_start_all_sat; r_start_sat];
   v_start_all_sat = [v_start_all_sat; v_start_sat];      
   
   for k = 1:size(r_start_all_sat,1)
       r_s = r_start_all_sat(k,:)
       sno = length(r_s)/3;
       
       for s = 1:sno
           sat_trx_pos_eci(s,1) = r_s(1,(3*s-2));
           sat_trx_pos_eci(s,2) = r_s(1,(3*s-1));
           sat_trx_pos_eci(s,3) = r_s(1,(3*s-0));
       end
   end   
   run('capacityVSdistance.m');
   
   xx(ss) = side/1000; % converting intersatellite distance meter to Km.
   y2(ss) = R;     
end

figure()
hold on;
plot3(sat_trx_pos_eci(:,1),sat_trx_pos_eci(:,2),sat_trx_pos_eci(:,3),'bv','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of transmitter antennas on 2 LEO satellites')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('ECI coordinates of Tx')
grid on;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% triangular configuration
r_start_all_sat = []; % Initialize to null.
v_start_all_sat = []; % Initialize to null.
for ss = 1:length(inter_satellite_distance)
   side = inter_satellite_distance(ss)  
   run('gco_hcw_triangle.m');
   run('init_start.m');
   r_start_all_sat = [r_start_all_sat; r_start_sat];
   v_start_all_sat = [v_start_all_sat; v_start_sat];      
   
   for k = 1:size(r_start_all_sat,1)
       r_s = r_start_all_sat(k,:)
       sno = length(r_s)/3;
       
       for s = 1:sno
           sat_trx_pos_eci(s,1) = r_s(1,(3*s-2));
           sat_trx_pos_eci(s,2) = r_s(1,(3*s-1));
           sat_trx_pos_eci(s,3) = r_s(1,(3*s-0));
       end
   end   
   run('capacityVSdistance.m');
   
   xx(ss) = side/1000; % converting intersatellite distance meter to Km.
   y3(ss) = R;  
   
end
figure()
hold on;
plot3(sat_trx_pos_eci(:,1),sat_trx_pos_eci(:,2),sat_trx_pos_eci(:,3),'bv','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of transmitter antennas on 3 LEO satellites')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('ECI coordinates of Tx')
grid on;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Square configuration
r_start_all_sat = []; % Initialize to null.
v_start_all_sat = []; % Initialize to null.
for ss = 1:length(inter_satellite_distance)
   side = inter_satellite_distance(ss)  
   run('gco_hcw_square.m');
   run('init_start.m');
   r_start_all_sat = [r_start_all_sat; r_start_sat];
   v_start_all_sat = [v_start_all_sat; v_start_sat];      
   
   for k = 1:size(r_start_all_sat,1)
       r_s = r_start_all_sat(k,:)
       sno = length(r_s)/3;
       
       for s = 1:sno
           sat_trx_pos_eci(s,1) = r_s(1,(3*s-2));
           sat_trx_pos_eci(s,2) = r_s(1,(3*s-1));
           sat_trx_pos_eci(s,3) = r_s(1,(3*s-0));
       end
   end   
   run('capacityVSdistance.m');
   
   xx(ss) = side/1000; % converting intersatellite distance meter to Km.
   y4(ss) = R;  
   
end
figure()
hold on;
plot3(sat_trx_pos_eci(:,1),sat_trx_pos_eci(:,2),sat_trx_pos_eci(:,3),'bv','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of transmitter antennas on 4 LEO satellites')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('ECI coordinates of Tx')
grid on;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Pentagon configuration
r_start_all_sat = []; % Initialize to null.
v_start_all_sat = []; % Initialize to null.
for ss = 1:length(inter_satellite_distance)
   side = inter_satellite_distance(ss)  
   run('gco_hcw_pentagon.m');
   run('init_start.m');
   r_start_all_sat = [r_start_all_sat; r_start_sat];
   v_start_all_sat = [v_start_all_sat; v_start_sat];      
   
   for k = 1:size(r_start_all_sat,1)
       r_s = r_start_all_sat(k,:)
       sno = length(r_s)/3;
       
       for s = 1:sno
           sat_trx_pos_eci(s,1) = r_s(1,(3*s-2));
           sat_trx_pos_eci(s,2) = r_s(1,(3*s-1));
           sat_trx_pos_eci(s,3) = r_s(1,(3*s-0));
       end
   end   
   run('capacityVSdistance.m');
   
   xx(ss) = side/1000; % converting intersatellite distance meter to Km.
   y5(ss) = R;     
end
figure()
hold on;
plot3(sat_trx_pos_eci(:,1),sat_trx_pos_eci(:,2),sat_trx_pos_eci(:,3),'bv','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of transmitter antennas on 5 LEO satellites')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('ECI coordinates of Tx')
grid on;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hexagon configuration
r_start_all_sat = []; % Initialize to null.
v_start_all_sat = []; % Initialize to null.
for ss = 1:length(inter_satellite_distance)
   side = inter_satellite_distance(ss)  
   run('gco_hcw_hexagon.m');
   run('init_start.m');
   r_start_all_sat = [r_start_all_sat; r_start_sat];
   v_start_all_sat = [v_start_all_sat; v_start_sat];      
   
   for k = 1:size(r_start_all_sat,1)
       r_s = r_start_all_sat(k,:)
       sno = length(r_s)/3;
       
       for s = 1:sno
           sat_trx_pos_eci(s,1) = r_s(1,(3*s-2));
           sat_trx_pos_eci(s,2) = r_s(1,(3*s-1));
           sat_trx_pos_eci(s,3) = r_s(1,(3*s-0));
       end
   end   
   run('capacityVSdistance.m');
   
   xx(ss) = side/1000; % converting intersatellite distance meter to Km.
   y6(ss) = R;     
end
figure()
hold on;
plot3(sat_trx_pos_eci(:,1),sat_trx_pos_eci(:,2),sat_trx_pos_eci(:,3),'bv','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of transmitter antennas on 6 LEO satellites')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('ECI coordinates of Tx')
grid on;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots the capacity vs inter satellite distance
figure();
hold on;

plot(xx,y2);
xlabel("Inter satellite distance (Km)");
ylabel("Achievable Channel Rate (R) (bps/Hz)");
title("Achievable Channel rate vs Inter-satellite distance");
grid on;
plot(xx,y3);
plot(xx,y4);
plot(xx,y5);
plot(xx,y6);
legend('Nt=1','Nt=3','Nt=4','Nt=5','Nt=6');


