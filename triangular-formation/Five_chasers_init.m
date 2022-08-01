clear all
close all
% initialize orbit calculation
% postion in metres, velocity in m/s
% SECOND satellite is the target/leader satellite
% FIRST, THIRD, FOURTH satellites are in triangular formation
r_start = [6.8705e+06 463.4268 -787.6687 ...
           6871000 0 0 ...
           6.8706e+06 -721.8017 -625.0976 ...
           6.8712e+06 -909.5208 401.3407 ...
           6.8715e+06 159.6820 873.1389 ...
           6.8711e+06 1.0082e+03 138.2919];
       
v_start = [-0.2568 7.6168e+03 0.4449 ...
           0 7616.269 0 ...
           0.4000 7.6167e+03  -0.6929 ...
           0.5041 7.6160e+03 -0.8731 ...
           -0.0885 7.6157e+03 0.1533 ...
           -0.5588 7.6162e+03 0.9678];
       
% start model
sim('model_plusacc2.slx')
 
% Plot orbit
plot(position(:,1),position(:,2))
hold on
plot(position(:,4),position(:,5))
plot(position(:,7),position(:,8))
plot(position(:,10),position(:,11))
plot(position(:,13),position(:,14))
plot(position(:,16),position(:,17))
rotate3d
title('orbit position ECI')
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal

% relative movement
% Pick a reference sat
sat_ref = 2;
anz = size(r_start);
sim_vec = size(position);
sim_points=sim_vec(1);
for i = 1 : anz(2)/3
sat_num=i;
for j = 1: sim_points
   delta_x(sat_num,j) = position(j,(sat_num-1)*3+1) - position(j,(sat_ref-1)*3+1); 
   delta_y(sat_num,j) = position(j,(sat_num-1)*3+2) - position(j,(sat_ref-1)*3+2); 
   delta_z(sat_num,j) = position(j,(sat_num-1)*3+3) - position(j,(sat_ref-1)*3+3); 
   % delta_x(j,i) = position(j,(i-1)*3+1) - position(j,1+(sat_ref-1)*3);
   % delta_y(j,i) = position(j,(i-1)*3+2) - position(j,1+(sat_ref-1)*3+1);
   % delta_z(j,i) = position(j,(i-1)*3+3) - position(j,1+(sat_ref-1)*3+2);
end
end
figure
hold on;
for i = 1 : anz(2)/3
plot3(delta_x(i,:),delta_y(i,:),delta_z(i,:));
title('relative position w.r.t. second satellite (leader) LVLH frame')
end
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')