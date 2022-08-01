%%%cartwheel2_GSfixedcoordinates.m%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts ECI of satellites to GS fixed coordinate system %

%side =inter_satellite_distance ; % meters, inter satellite distance is 2Km
run('gco_hcw_cartwheel.m');
run('init_start.m');
Nt = length(r_start_sat)/3 ; % number of transmitting antennas in space

ant_loc(1,1) = ant_location(1,1)
ant_loc(1,2) = ant_location(1,2)
ant_loc(1,3) = ant_location(1,3)

%%%% translation of axes %%%%
r_start_sat_gs=[];
for e = 1:(length(r_start_sat)/3)
    r_start_sat_gs(1,(3*e-2)) = r_start_sat(1,(3*e-2)) - ant_loc(1,1);
    r_start_sat_gs(1,(3*e-1)) = r_start_sat(1,(3*e-1)) - ant_loc(1,2);
    r_start_sat_gs(1,(3*e-0)) = r_start_sat(1,(3*e-0)) - ant_loc(1,3);    
end