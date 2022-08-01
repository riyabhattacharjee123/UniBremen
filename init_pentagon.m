close all
% initialize orbit calculation
% postion in metres, velocity in m/s
% SECOND satellite is the target/leader satellite
% FIRST, THIRD, FOURTH satellites are in triangular formation

% Data of chief
r_0 = [6871000 0 0];
v_0 = [0 sqrt(mu/r_0(1)) 0];


% Data of deputies
aux_p = [x_c y_c z_c].';
aux_v = [xdot_c ydot_c zdot_c].';
aux_p = aux_p(:).';
aux_v = aux_v(:).';

% Data of all satellites
r_start = [aux_p(1:3) r_0 aux_p(4:end)]
v_start = [aux_v(1:3) v_0 aux_v(4:end)]

%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Simulation of ideal satellite orbits
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
       
% start model
sim('model.slx')
 
% Plot orbit
figure();
hold on;
plot(position(:,1),position(:,2))
hold on;
plot(position(:,4),position(:,5))
plot(position(:,7),position(:,8))
plot(position(:,10),position(:,11))
plot(position(:,13),position(:,14))
plot(position(:,16),position(:,17))
rotate3d
title('orbit position ECI')
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
legend('chaser-1','leader','chaser-2','chaser-3','chaser-4','chaser-5')

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
axis equal
legend('chaser-1','chaser-2','chaser-3','chaser-4','chaser-5')


%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Simulations of orbits with added Perturbation Accelerations
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

% start model
sim('model_plusacc4.slx')

% Plot orbit
figure();
hold on;
plot(position(:,1),position(:,2))
hold on;
plot(position(:,4),position(:,5))
plot(position(:,7),position(:,8))
plot(position(:,10),position(:,11))
plot(position(:,13),position(:,14))
plot(position(:,16),position(:,17))
rotate3d
title('orbit position ECI (Adding Perturbations)')
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
legend('chaser-1','leader','chaser-2','chaser-3','chaser-4','chaser-5')

% relative movement
%  reference sat 2

for i = 1 : anz(2)/3
sat_num=i;
for j = 1: sim_points
   delta_x_a(sat_num,j) = position(j,(sat_num-1)*3+1) - position(j,(sat_ref-1)*3+1); 
   delta_y_a(sat_num,j) = position(j,(sat_num-1)*3+2) - position(j,(sat_ref-1)*3+2); 
   delta_z_a(sat_num,j) = position(j,(sat_num-1)*3+3) - position(j,(sat_ref-1)*3+3);    
end
end
figure();
hold on;
for i = 1 : anz(2)/3
plot3(delta_x_a(i,:),delta_y_a(i,:),delta_z_a(i,:));
title('relative position w.r.t. second satellite (leader) LVLH frame (Adding Perturbations)')
end
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')
axis equal
legend('chaser-1','leader','chaser-2','chaser-3','chaser-4','chaser-5')
