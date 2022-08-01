% Convert geodetic latitude, longitude, altitude (LLA) coordinates to 
% Earth-centered inertial (ECI) coordinates

% Antenna locations
run('receiver_antenna_locations.m')

% Convert ECI (cartesian) into Polar coordinates for satellite

% Satellite position taken from triangle simulation Satellite position
%r_start = [6871288.67513460,1000.00000000000,500.000000000000,6871000,0,0,6871288.67513460,-1000.00000000000,500.000000000000,6870422.64973081,-2.12115047744981e-13,-1000];

for ri = 1:length(r_start)
    disp(r_start)
    snum = length(r_start)/3;
    for s = 1:snum        
        if s ~= 2
            trx_pos_eci(s,1) = r_start(1,(3*s-2));
            trx_pos_eci(s,2) = r_start(1,(3*s-1));
            trx_pos_eci(s,3) = r_start(1,(3*s-0));
        end
    end    
end
trx_pos_eci([2],:) = [];
disp(trx_pos_eci);

% convert Transmiter position ECI to Polar coordinates
trx_pos_pol = []; % [r theta phi ; r theta phi];
for r1=1:size(trx_pos_eci,1)
    %fprintf("This is row:: %d \n",r1 )
     %fprintf("This are row elements : %f \n",trx_pos_eci(r1,:) )
     x_sq = trx_pos_eci(r1,1)^2;
      y_sq = trx_pos_eci(r1,2)^2;
       z_sq = trx_pos_eci(r1,3)^2;
       
     r = sqrt(x_sq+y_sq+z_sq); % meters radius in polar coordinates   
     % theta polar cood.
     theta = atan((sqrt(x_sq+y_sq))/trx_pos_eci(r1,3)); %radians 
     % phi polar cood.
     phi = atan(trx_pos_eci(r1,2)/trx_pos_eci(r1,1)) ; % radians
     
     % assign polar cood to matrix
     trx_pos_pol(r1,1) = r;  % meters radius in polar coordinates   
     trx_pos_pol(r1,2) = theta ; % radians
     trx_pos_pol(r1,3) = phi ;   % radians  
    
end

%disp(trx_pos_pol)

% convert Receiver position ECI to Polar coordinates
radius_earth = 6371000 ; %metres
rcx_pos_pol = [radius_earth pi/2 0]; % [radius_earth theta phi];     
     % assign polar cood to matrix
rcx_pos_pol(1) = radius_earth;  % meters radius in polar coordinates   
rcx_pos_pol(2) = atan((sqrt(x_sq+y_sq))/trx_pos_eci(r1,3)); %radians 
rcx_pos_pol(3) =  atan(trx_pos_eci(r1,2)/trx_pos_eci(r1,1)) ; % radians 

%disp(rcx_pos_pol)


