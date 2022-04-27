% Convert geodetic latitude, longitude, altitude (LLA) coordinates to 
% Earth-centered inertial (ECI) coordinates

pi = 22/7;
% Number of Earth station receivers
Nr=3;

% Receiver Antenna locations in Lat., Lon., Alt.
rcvr_1=[53.073635 8.806422 18]; % Bremen Uni 53°N , 8.8°E, 18m 
rcvr_2=[53.139403 8.677612 17]; % Bremen Industriehäfen
rcvr_3=[53.046563 8.744246 15]; % Bremen Roland Center

rcvr_pos_lla = [rcvr_1;rcvr_2;rcvr_3];

rcvr_pos_eci_1 = lla2eci(rcvr_1,utc);
rcvr_pos_eci_2 = lla2eci(rcvr_2,utc);
rcvr_pos_eci_3 = lla2eci(rcvr_3,utc);

rcvr_pos_eci = [rcvr_pos_eci_1;
    rcvr_pos_eci_2;
    rcvr_pos_eci_3];


% Convert ECI (cartesian) into Polar coordinates for satellite

% Satellite position taken from triangle simulation Satellite position
trx_pos_eci = [6.741288675134595e+06 1.0e+03 4.999999999999999e+02;
    6.741288675134595e+06 -1.000000000000000e+03 4.999999999999999e+02;
    6.740422649730810e+06 -2.121150477449814e-13 -1000];

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


