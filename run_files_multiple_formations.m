clear;
close all;
clear;
clc;
run('User_input_data.m');
run('receiver_antenna_locations.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 SATELLITES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LEADER FOLLOWER FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('LF2_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('sat_2_steeringVectors.m');
    run('channel_rate3Dmax.m');    
    
    Y2LF(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PENDULUM FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('pendulum2_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('sat_2_steeringVectors.m');
    run('channel_rate3Dmax.m');    
    
    Y2Pen(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CARTWHEEL FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('cartwheel2_GSfixedcoordinates.m');
    run('cartesianTOspherical_r_start_sat.m'); 
    run('sat_2_steeringVectors.m');
    run('channel_rate3Dmax.m');    
    
    Y2Cart(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3 SATELLITES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LEADER FOLLOWER FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('LF3_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors.m');
    run('channel_rate3Dmax.m');     
    
    Y3LF(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PENDULUM FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('pendulum3_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors.m'); % 3 satellites
    run('channel_rate3Dmax.m');    
    
    Y3Pen(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CARTWHEEL FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('cartwheel3_GSfixedcoordinates.m');
    run('cartesianTOspherical_r_start_sat.m'); 
    run('steeringVectors.m'); % 3 satellites
    run('channel_rate3Dmax.m');     
    
    Y3Cart(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRIANGLE FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('GSfixedcoordinates.m');
    run('cartesianTOspherical_r_start_sat.m'); 
    run('steeringVectors.m');
    run('channel_rate3Dmax.m');        
    Y3Tri(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4 SATELLITES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LEADER FOLLOWER FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('LF4_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors_4.m');
    run('channel_rate3Dmax.m');     
    
    Y4LF(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRIANGLE FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('Tri_4_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors_4.m');
    run('channel_rate3Dmax.m');     
    
    Y4Tri(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SQUARE FORMATION %%%
for isd = 1: length(inter_satellite_distance)
     side = inter_satellite_distance(isd)
    run('GSfixedcoordinates_4.m');
    run('cartesianTOspherical_r_start_sat.m'); 
    run('steeringVectors_4.m');
    run('channel_rate3Dmax.m');     
    
    Y4Sq(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5 SATELLITES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SQUARE FORMATION %%%
for isd = 1: length(inter_satellite_distance)
    side = inter_satellite_distance(isd)
    run('Sq_5_GSfixedcoordinates.m');
    run('cartesianTOspherical.m'); 
    run('steeringVectors_5.m');
    run('channel_rate3Dmax.m');    
    
    Y5Sq(isd)= real(R_opt) ; 
    X1(isd)= inter_satellite_distance(isd)/1000;         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot graph
figure();
plot(X1,Y2LF,'--o');
hold on;
plot(X1,Y2Pen, 'r--o');
plot(X1,Y2Cart, 'b--o');
plot(X1,Y3LF);
plot(X1,Y3Pen);
plot(X1,Y3Cart);
plot(X1,Y3Tri);
plot(X1,Y4LF,'--^');
plot(X1,Y4Tri,'--^');
plot(X1,Y4Sq,'--^');
plot(X1,Y5Sq,'-*');

ylabel("Optimum Channel Rate (bps/Hz)");
xlabel("Inter satellite distance (Km)");
title("Distance between satellites vs R optimum");
legend('2 sat. LF','2 sat. Pendulum',...
    '2 sat. Cartwheel',...
    '3 sat. LF',...
    '3 sat. Pendulum',...
    '3 sat. Cartwheel',...
    '3 sat. Triangle',...
    '4 sat. LF',...
    '4 sat. Triangle',...
    '4 sat. Square',...
    '5 sat. square')
grid on;