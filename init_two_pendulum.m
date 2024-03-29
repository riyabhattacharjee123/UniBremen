close all
% initialize orbit calculation
% postion in metres, velocity in m/s
% SECOND satellite is the target/leader satellite
% FIRST satellite is the trailling/chaser satellite
r_0 = [6741000 0 0] ;
v_0 = [0 sqrt(mu/r_0(1)) 0];

% Data of chaser
aux_p = [x_c y_c z_c].';
aux_v = [xdot_c ydot_c zdot_c].';
aux_p = aux_p(:).';
aux_v = aux_v(:).';

% Data of all satellites
r_start = [aux_p(1:3) r_0 aux_p(4:end)]
v_start = [aux_v(1:3) v_0 aux_v(4:end)]

% start model
sim('model_plusacc.slx')
 
% Plot orbit
plot(position(:,1),position(:,2))
hold on
plot(position(:,4),position(:,5))
%plot(position(:,7),position(:,8))
%plot(position(:,10),position(:,11))
%plot(position(:,13),position(:,14))
rotate3d
title('orbit position in ECI frame')
xlabel('ECI-x'), ylabel('ECI-y'), zlabel('ECI-z')
axis equal
legend('chaser-1','leader')

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
title('relative position w.r.t. second satellite (leader) in LVLH frame')
end
rotate3d
xlabel('LVLH-x'), ylabel('LVLH-y'), zlabel('LVLH-z')
legend('chaser orbit in pendulum formation')