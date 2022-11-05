%%get_eci_loc_sat.m%%
% gets the ECI location of Satellite when the Receiver ECI is given

% 01/17/2022 10:20:36 UTC (MM/DD/YYYY HH:MM:SS)
utc = [2022 1 17 10 20 36];
% Receiver Antenna locations in Lat., Lon., Alt.
rcvr_1=[53.10657108482692, 8.850392209543823 18]; % Bremen Uni 53°N , 8.8°E, 18m 
Nrx= 8;
ant_location = lla2eci(rcvr_1,utc);
rcvr_pos_eci = ant_location;
point1=rcvr_pos_eci(1,:);

[gs_azimuth,gs_elevation,gs_r] = cart2sph(rcvr_pos_eci(1,1),rcvr_pos_eci(1,2),rcvr_pos_eci(1,3));

sat_eci(1,1) = rcvr_pos_eci(1,1) + d_r*cosd(theta_el_deg)*cosd(phi_az_deg);
sat_eci(1,2) = rcvr_pos_eci(1,2) + d_r*sind(theta_el_deg);
sat_eci(1,3) = rcvr_pos_eci(1,3) + d_r*cosd(theta_el_deg)*sind(360-phi_az_deg);


