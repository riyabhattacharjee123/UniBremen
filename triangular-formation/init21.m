clear all
close all
% initialize orbit calculation
% postion in metres, velocity in m/s
% SECOND satellite is the target/leader satellite
% FIRST, THIRD, FOURTH satellites are in triangular formation
r_start = [6.8713e+06 999.9999 500.0002 ...
           6871000 0 0 ...
           6.8713e+06 -999.9978 500.0034 ...
           6.8704e+06 0.0012 -1.0000e+03];
       
v_start = [-0.5542 7.6159e+03 0.9600 ...
           0 7616.269 0 ...
            0.5542 7.6159e+03  -0.9600 ...
           -6.5253e-07 7.6169e+03 1.1302e-06];
       
% start model
sim('model_plusacc2.slx')
 
% Plot orbit
plot(position(:,1),position(:,2))
hold on
plot(position(:,4),position(:,5))
plot(position(:,7),position(:,8))
plot(position(:,10),position(:,11))
%plot(position(:,13),position(:,14))
rotate3d
title('orbit position')
xlabel('x'), ylabel('y'), zlabel('z')
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
title('relative position w.r.t. second satellite (leader)')
end
rotate3d
xlabel('x'), ylabel('y'), zlabel('z')