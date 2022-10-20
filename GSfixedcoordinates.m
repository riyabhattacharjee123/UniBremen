%%%GSfixedcoordinates.m%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts ECI of satellites to GS fixed coordinate system %

%side =inter_satellite_distance ; % meters, inter satellite distance is 2Km
run('gco_hcw_triangle.m');
run('init_start.m');
Nt = length(r_start_sat)/3 ; % number of transmitting antennas in space
Ns = Nt ; % number of satellites in space

%%%% translation of axes %%%%
r_start_sat_gs=[];
ant_loc(1,1) = ant_location(1,1)
ant_loc(1,2) = ant_location(1,2)
ant_loc(1,3) = ant_location(1,3)

for e = 1:(length(r_start_sat)/3)
    % GS centered x-axis
    r_start_sat_gs(1,(3*e-2)) = r_start_sat(1,(3*e-2)) - ant_loc(1,1); 
    
    % GS centered y-axis
    r_start_sat_gs(1,(3*e-1)) = r_start_sat(1,(3*e-1)) - ant_loc(1,2);
    
    % GS centered z-axis
    r_start_sat_gs(1,(3*e-0)) = r_start_sat(1,(3*e-0)) - ant_loc(1,3);    
end