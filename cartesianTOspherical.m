% cartesianTOspherical.m%
% converts cartesian to spherical coordinates %
%%%% translation of axes %%%%
% convert the GS centered cartesian coordinates to GS centered spherical
% cood.

r_start_sat_gs_spherical=[]; % meters, rad, rad. Distance,
%polar angle theta, azimuth angle phi

for f = 1:(length(r_start)/3)
    
    x_dash =  r_start_sat_gs(1,(3*f-2));
    z_dash = r_start_sat_gs(1,(3*f-1));
    y_dash = r_start_sat_gs(1,(3*f-0));
    
    % radius from GS to satellite distance
    r_start_sat_gs_spherical(1,(3*f-2))= sqrt(x_dash^2 ...
             + y_dash^2 ...
             + z_dash^2 ); % meters
    
    %polar/inclination angle theta from GS centered y-axis
    r_start_sat_gs_spherical(1,(3*f-1)) = atan( sqrt(x_dash^2+z_dash^2)...
                                          / y_dash ); % radians
                                      
                                      
    
                                      
    %r_start_sat_gs_spherical(1,(3*f-2)) = pi/2 ...
       %                                   - r_start_sat_gs_spherical...
      %                                    (1,(3*f-3)) ; % radians
                                      %elevation angle from x-z plane
                                      
    
    % azimuth angle from GS centered x-axis on xz-plane                                  
    r_start_sat_gs_spherical(1,(3*f-0)) = atan(z_dash/x_dash); % radians
        
    %r_start_sat_gs_spherical(1,(3*f-0)) = 2*pi - atan(z_dash/x_dash); %radians,
    %azimuth angle phi taken from x-axis till z-axis
end

