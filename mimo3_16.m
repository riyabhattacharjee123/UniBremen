% Antenna locations
run('receiver_antenna_locations.m')

% satellite positions
inter_satellite_distance = (50:100:50000); % meters

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

figure();
hold on;

plot(xx,y3);
xlabel("Inter satellite distance (Km)");
ylabel("Achievable Channel Rate (R) (bps/Hz)");
title("Achievable Channel rate vs Inter-satellite distance");
grid on;
legend('3X16 MIMO')