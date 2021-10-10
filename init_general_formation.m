close all
% initialize orbit calculation
% postion in metres, velocity in m/s
% SECOND satellite is the target/leader satellite
% FIRST, THIRD, FOURTH satellites are in triangular formation

%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Run cwgco_code.m to get the data of the chief and deputy satellites
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

% Data of chief
r_0 = [x_leader y_leader z_leader];
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

% Second satellite is the leader satellite
sat_ref = 2;
% calculate the size of input vector
anz = size(r_start);
% calculate the simulation vector
sim_vec = size(position);
% calculat the simulation points
sim_points=sim_vec(1);
 
% Plot absolute orbit in ECI
% ideal scenario, without external perturbations
figure();
hold on;
N = (length(r_start)/3)
for c = 1 : N
plot3(position(:,(3*c-2)),position(:,(3*c-1)),position(:,(3*c-0)));
hold on;
title('orbit position ECI');
end
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
Legend1=cell(N,1);
for iter=1:N
    if iter==2
        Legend1{iter}=strcat('Leader Satellite- ', num2str(iter));
    else
        Legend1{iter}=strcat('Deputy satellite- ', num2str(iter));
    end
end
legend(Legend1)
%legend('chaser-1','leader','chaser-2','chaser-3')

% relative movement with respect to the leader satellite
% ideal scenario, without external perturbations
for i = 1 : anz(2)/3
sat_num=i;
for j = 1: sim_points
   delta_x(sat_num,j) = position(j,(sat_num-1)*3+1) - position(j,(sat_ref-1)*3+1); 
   delta_y(sat_num,j) = position(j,(sat_num-1)*3+2) - position(j,(sat_ref-1)*3+2); 
   delta_z(sat_num,j) = position(j,(sat_num-1)*3+3) - position(j,(sat_ref-1)*3+3); 
end
end
figure();
hold on;
for i = 1 : anz(2)/3
plot3(delta_x(i,:),delta_y(i,:),delta_z(i,:));
title('relative position w.r.t. second satellite (leader) LVLH frame')
end
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')
axis equal
Legend2=cell((N-1),1);
for iter=1:(N-1)
    Legend2{iter}=strcat('Deputy satellite- ', num2str(iter));
end
legend(Legend2)
%legend('chaser-1','chaser-2','chaser-3')

%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Simulations of orbits with added Perturbation Accelerations
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

% start model
% ( Perturbations added )
sim('model_plusacc4.slx')

% Plot absolute orbit in ECI
% ( Perturbations added )
figure();
hold on;
for c = 1 : (length(r_start)/3)
plot3(position_perturbed(:,(3*c-2)),position_perturbed(:,(3*c-1)),position_perturbed(:,(3*c-0)));
hold on;
title('orbit position ECI (Adding Perturbations)')
end
rotate3d;
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
Legend3=cell(N,1);
for iter=1:N
    if iter==2
        Legend3{iter}=strcat('Leader Satellite- ', num2str(iter));
    else
        Legend3{iter}=strcat('Deputy satellite- ', num2str(iter));
    end
end
legend(Legend3)
%legend('chaser-1','leader','chaser-2','chaser-3')

% relative movement with respect to the leader satellite
% ( Perturbations added )
% Pick a reference sat
sat_ref1 = 2;
anz = size(r_start);
sim_vec1 = size(position_perturbed);
sim_points1=sim_vec1(1);
for i = 1 : anz(2)/3
sat_num1=i;
for j = 1: sim_points1
   delta_x_a(sat_num1,j) = position_perturbed(j,(sat_num1-1)*3+1) - position_perturbed(j,(sat_ref1-1)*3+1); 
   delta_y_a(sat_num1,j) = position_perturbed(j,(sat_num1-1)*3+2) - position_perturbed(j,(sat_ref1-1)*3+2); 
   delta_z_a(sat_num1,j) = position_perturbed(j,(sat_num1-1)*3+3) - position_perturbed(j,(sat_ref1-1)*3+3);    
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
Legend4=cell((N-1),1);
for iter=1:(N-1)
    Legend4{iter}=strcat('Deputy satellite- ', num2str(iter));
end
legend(Legend4)
%legend('chaser-1','chaser-2','chaser-3')


%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Plotting the simulation differences over time
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

