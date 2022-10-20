%close all
% init_general_formation.m %
% initialize orbit calculation
% postion in metres, velocity in m/s, acceleration in m/s^2
% SECOND satellite is the target/leader satellite
% Other satellites in position 1, 3 ,4 etc (as applicable) are Deputies

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
sim('model_ideal_2017b.slx')% this refers to the model without perturbations

% Second satellite is the leader satellite
sat_ref = 2;
% calculate the size of input vector
anz = size(r_start);
% calculate the simulation vector
sim_vec = size(position);
% calculate the simulation points
sim_points=sim_vec(1);
 
% Plot absolute orbit in ECI
% ideal scenario, without external perturbations
figure();
hold on;
N = (length(r_start)/3)
for c = 1 : N
plot3(position(:,(3*c-2)),position(:,(3*c-1)),position(:,(3*c-0)));
set(gca,'FontSize',15);
hold on;
title('Orbit position ECI');
end
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
Legend1=cell(N,1);
for iter=1:N
    if iter==2
        Legend1{iter}=strcat('Leader/reference- ', num2str(iter));
    else
        Legend1{iter}=strcat('Deputy- ', num2str(iter));
    end
end
legend(Legend1)

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
set(gca,'FontSize',15);
title('Relative position w.r.t. leader satellite LVLH frame')
end
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')
axis equal
Legend2=cell((N-1),1);
for iter=1:(N-1)
    Legend2{iter}=strcat('Deputy- ', num2str(iter));
end
legend(Legend2)

%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Simulations of orbits with added Perturbation Accelerations
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

% start model
% ( Perturbations added )
sim('model_plusacc_2017b.slx') % this refers to the model with added perturbations

% Plot absolute orbit in ECI
% ( Perturbations added )
figure();
hold on;
for c = 1 : (length(r_start)/3)
plot3(position_perturbed(:,(3*c-2)),position_perturbed(:,(3*c-1)),position_perturbed(:,(3*c-0)));
set(gca,'FontSize',15);
hold on;
title('Orbit position ECI (Perturbed)')
end
rotate3d;
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
axis equal
Legend3=cell(N,1);
for iter=1:N
    if iter==2
        Legend3{iter}=strcat('Leader/reference- ', num2str(iter));
    else
        Legend3{iter}=strcat('Deputy- ', num2str(iter));
    end
end
legend(Legend3)

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
set(gca,'FontSize',15);
title('Relative position w.r.t. leader satellite LVLH frame (Perturbed)')
end
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')
axis equal
Legend4=cell((N-1),1);
for iter=1:(N-1)
    Legend4{iter}=strcat('Deputy- ', num2str(iter));
end
legend(Legend4)

%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%
% Plotting the simulation differences over time
% Plotting timeseries in ECI-x , ECI-y, ECI-z direction
% Plotting the differences in the ideal and perturbed paths
% position and velocity both considered
%------------------------------------------------------------------------------------------------%
%------------------------------------------------------------------------------------------------%

% Position in ECI-x , ECI-y, ECI-z direction

for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x position
    nexttile
    pos_x_1 = position(:,(3*s-2));
    ts_pos_x_1 = timeseries(pos_x_1,1:sim_points);
    ts_pos_x_1.Time = ts_pos_x_1.Time - ts_pos_x_1.Time(1); 
    plot(ts_pos_x_1)
    set(gca,'FontSize',15);
    hold on;
    pos_x_2 = position_perturbed(:,(3*s-2));
    ts_pos_x_2 = timeseries(pos_x_2,1:sim_points1);
    ts_pos_x_2.Time = ts_pos_x_2.Time - ts_pos_x_2.Time(1); 
    plot(ts_pos_x_2)
    set(gca,'FontSize',15);
    
    pos_diff_x = pos_x_1 - pos_x_2;
    ts_pos_diff_x =  timeseries(pos_diff_x,1:sim_points);
    ts_pos_diff_x.Time = ts_pos_diff_x.Time - ts_pos_diff_x.Time(1); 
    plot(ts_pos_diff_x)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Position in ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (metres)')
    legend('Ideal path','Perturbed path','position difference');

    % ECI-y position
    nexttile
    pos_y_1 = position(:,(3*s-1));
    ts_pos_y_1 = timeseries(pos_y_1,1:sim_points);
    ts_pos_y_1.Time = ts_pos_y_1.Time - ts_pos_y_1.Time(1); 
    plot(ts_pos_y_1)
    set(gca,'FontSize',15);
    hold on;
    pos_y_2 = position_perturbed(:,(3*s-1));
    ts_pos_y_2 = timeseries(pos_y_2,1:sim_points1);
    ts_pos_y_2.Time = ts_pos_y_2.Time - ts_pos_y_2.Time(1); 
    plot(ts_pos_y_2)
    set(gca,'FontSize',15);
    
    pos_diff_y = pos_y_1 - pos_y_2;
    ts_pos_diff_y =  timeseries(pos_diff_y,1:sim_points);
    ts_pos_diff_y.Time = ts_pos_diff_y.Time - ts_pos_diff_y.Time(1); 
    plot(ts_pos_diff_y)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Position in ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (metres)')
    legend('Ideal path','Perturbed path','position difference');
    
    % ECI-z position
    nexttile
    pos_z_1 = position(:,(3*s-0));
    ts_pos_z_1 = timeseries(pos_z_1,1:sim_points);
    ts_pos_z_1.Time = ts_pos_z_1.Time - ts_pos_z_1.Time(1); 
    plot(ts_pos_z_1)
    set(gca,'FontSize',15);
    hold on;
    pos_z_2 = position_perturbed(:,(3*s-0));
    ts_pos_z_2 = timeseries(pos_z_2,1:sim_points1);
    ts_pos_z_2.Time = ts_pos_z_2.Time - ts_pos_z_2.Time(1); 
    plot(ts_pos_z_2)
    set(gca,'FontSize',15);
    
    pos_diff_z = pos_z_1 - pos_z_2;
    ts_pos_diff_z =  timeseries(pos_diff_z,1:sim_points);
    ts_pos_diff_z.Time = ts_pos_diff_z.Time - ts_pos_diff_z.Time(1); 
    plot(ts_pos_diff_z)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Position in ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (metres)')
    legend('Ideal path','Perturbed path','position difference');
    
    % only differences in position
    figure();
    
    nexttile
    plot(ts_pos_diff_x)    
    set(gca,'FontSize',15);
    title('Position difference: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (metres)')
    legend('position difference');
    hold on;
    
    nexttile
    plot(ts_pos_diff_y)  
    set(gca,'FontSize',15);
    title('Position difference: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (metres)')
    legend('position difference');
        
    nexttile
    plot(ts_pos_diff_z)   
    set(gca,'FontSize',15);
    title('Position difference: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (metres)')
    legend('position difference');     
    
    nexttile
     % Magnitude of Unperturbed Position
     position_mag =  sqrt(pos_x_1.^2+pos_y_1.^2+pos_z_1.^2);
     plot(position_mag);
     set(gca,'FontSize',15);
     hold on;
      % Magnitude of Perturbed Position
     position_pert_mag = sqrt(pos_x_2.^2+pos_y_2.^2+pos_z_1.^2);
     plot(position_pert_mag);
     set(gca,'FontSize',15);
     title('Magnitude of the Position Vector for Satellite-',num2str(s));
     xlabel('Time(seconds)'), ylabel('Position (metres)');
     legend('Ideal path magnitude','Perturbed path magnitude');
     xlim([0 10000]);   
     
end


% Velocity in ECI-x, ECI-y, ECI-z direction

for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Velocity
    nexttile
    vel_x_1 = velocity(:,(3*s-2));
    ts_vel_x_1 = timeseries(vel_x_1,1:sim_points);
    ts_vel_x_1.Time = ts_vel_x_1.Time - ts_vel_x_1.Time(1); 
    plot(ts_vel_x_1)
    set(gca,'FontSize',15);
    hold on;
    vel_x_2 = velocity_perturbed(:,(3*s-2));
    ts_vel_x_2 = timeseries(vel_x_2,1:sim_points1);
    ts_vel_x_2.Time = ts_vel_x_2.Time - ts_vel_x_2.Time(1); 
    plot(ts_vel_x_2)
    set(gca,'FontSize',15);
    
    vel_diff_x = vel_x_1 - vel_x_2;
    ts_vel_diff_x =  timeseries(vel_diff_x,1:sim_points);
    ts_vel_diff_x.Time = ts_vel_diff_x.Time - ts_vel_diff_x.Time(1); 
    plot(ts_vel_diff_x)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Velocity in ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s)')
    legend('Ideal path','Perturbed path','velocity difference');

    % ECI-y velocity
    nexttile
    vel_y_1 = velocity(:,(3*s-1));
    ts_vel_y_1 = timeseries(vel_y_1,1:sim_points);
    ts_vel_y_1.Time = ts_vel_y_1.Time - ts_vel_y_1.Time(1); 
    plot(ts_vel_y_1)
    set(gca,'FontSize',15);
    hold on;
    vel_y_2 = velocity_perturbed(:,(3*s-1));
    ts_vel_y_2 = timeseries(vel_y_2,1:sim_points1);
    ts_vel_y_2.Time = ts_vel_y_2.Time - ts_vel_y_2.Time(1); 
    plot(ts_vel_y_2)
    set(gca,'FontSize',15);
    
    vel_diff_y = vel_y_1 - vel_y_2;
    ts_vel_diff_y =  timeseries(vel_diff_y,1:sim_points);
    ts_vel_diff_y.Time = ts_vel_diff_y.Time - ts_vel_diff_y.Time(1); 
    plot(ts_vel_diff_y)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Velocity in ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s)')
    legend('Ideal path','Perturbed path','velocity difference');
    
     % ECI-z velocity
    nexttile
    vel_z_1 = velocity(:,(3*s-0));
    ts_vel_z_1 = timeseries(vel_z_1,1:sim_points);
    ts_vel_z_1.Time = ts_vel_z_1.Time - ts_vel_z_1.Time(1); 
    plot(ts_vel_z_1)
    set(gca,'FontSize',15);
    hold on;
    vel_z_2 = velocity_perturbed(:,(3*s-0));
    ts_vel_z_2 = timeseries(vel_z_2,1:sim_points1);
    ts_vel_z_2.Time = ts_vel_z_2.Time - ts_vel_z_2.Time(1); 
    plot(ts_vel_z_2)
    set(gca,'FontSize',15);
    
    vel_diff_z = vel_z_1 - vel_z_2;
    ts_vel_diff_z =  timeseries(vel_diff_z,1:sim_points);
    ts_vel_diff_z.Time = ts_vel_diff_z.Time - ts_vel_diff_z.Time(1); 
    plot(ts_vel_diff_z)
    set(gca,'FontSize',15);
    
    title('Difference in Ideal and Perturbed orbit Velocity in ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s)')
    legend('Ideal path','Perturbed path','velocity difference');
    
    % only differences in velocity
    figure();
    
    nexttile
    plot(ts_vel_diff_x)   
    set(gca,'FontSize',15);
    title('Velocity difference: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s)')
    legend('velocity difference');
    hold on;
    
    nexttile
    plot(ts_vel_diff_y)    
    set(gca,'FontSize',15);
    title('Velocity difference: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s)')
    legend('velocity difference');
        
    nexttile
    plot(ts_vel_diff_z)    
    set(gca,'FontSize',15);
    title('Velocity difference: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s)')
    legend('velocity difference'); 
    
    nexttile
     % Magnitude of Unperturbed Velocity
     vel_mag =  sqrt(vel_x_1.^2+vel_y_1.^2+vel_z_1.^2);
     plot(vel_mag);
     set(gca,'FontSize',15);
     hold on;
      % Magnitude of Perturbed Velocity
     vel_pert_mag = sqrt(vel_x_2.^2+vel_y_2.^2+vel_z_1.^2);
     plot(vel_pert_mag);
     set(gca,'FontSize',15);
     title('Magnitude of the Velocity Vector for Satellite-',num2str(s));
     xlabel('Time(seconds)'), ylabel('Velocity (metres/second)');
     legend('Ideal velocity magnitude','Perturbed velocity magnitude');
     xlim([0 10000])
    
     
end

% Acceleration in ECI-x, ECI-y, ECI-z direction

for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)
    set(gca,'FontSize',15);
    hold on;
    acc_x_2 = acc_perturbed(:,(3*s-2));
    ts_acc_x_2 = timeseries(acc_x_2,1:sim_points1);
    ts_acc_x_2.Time = ts_acc_x_2.Time - ts_acc_x_2.Time(1); 
    plot(ts_acc_x_2)
    set(gca,'FontSize',15);
    
    acc_diff_x = acc_x_1 - acc_x_2;
    ts_acc_diff_x =  timeseries(acc_diff_x,1:sim_points);
    ts_acc_diff_x.Time = ts_acc_diff_x.Time - ts_acc_diff_x.Time(1); 
    plot(ts_acc_diff_x)
    set(gca,'FontSize',15);
    
    title('Ideal and Perturbed orbit Acceleration in ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Ideal path','Perturbed path','acceleration difference');

    % ECI-y Acceleration
    nexttile
    acc_y_1 = acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)
    set(gca,'FontSize',15);
    hold on;
    acc_y_2 = acc_perturbed(:,(3*s-1));
    ts_acc_y_2 = timeseries(acc_y_2,1:sim_points1);
    ts_acc_y_2.Time = ts_acc_y_2.Time - ts_acc_y_2.Time(1); 
    plot(ts_acc_y_2)
    
    acc_diff_y = acc_y_1 - acc_y_2;
    ts_acc_diff_y =  timeseries(acc_diff_y,1:sim_points);
    ts_acc_diff_y.Time = ts_acc_diff_y.Time - ts_acc_diff_y.Time(1); 
    plot(ts_vel_diff_y)
    set(gca,'FontSize',15);
    
    title('Ideal and Perturbed orbit Acceleration in ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Ideal path','Perturbed path','acceleration difference');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
    set(gca,'FontSize',15);
    hold on;
    acc_z_2 = acc_perturbed(:,(3*s-0));
    ts_acc_z_2 = timeseries(acc_z_2,1:sim_points1);
    ts_acc_z_2.Time = ts_acc_z_2.Time - ts_acc_z_2.Time(1); 
    plot(ts_acc_z_2)
    set(gca,'FontSize',15);
    
    acc_diff_z = acc_z_1 - acc_z_2;
    ts_acc_diff_z =  timeseries(acc_diff_z,1:sim_points);
    ts_acc_diff_z.Time = ts_acc_diff_z.Time - ts_acc_diff_z.Time(1); 
    plot(ts_acc_diff_z)
    set(gca,'FontSize',15);
    
    title('Ideal and Perturbed orbit Acceleration in ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Ideal path','Perturbed path','acceleration difference');
    
    % only differences in acceleration
    figure();
    
    nexttile
    plot(ts_acc_diff_x)    
    set(gca,'FontSize',15);
    title('Acceleration difference: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('acceleration difference');
    hold on;
    
    nexttile
    plot(ts_acc_diff_y)  
    set(gca,'FontSize',15);
    title('Acceleration difference: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('acceleration difference');
    
    nexttile
    plot(ts_acc_diff_z)  ;
    set(gca,'FontSize',15);
    title('Acceleration difference: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('acceleration difference');
    
    nexttile
     % Magnitude of Unperturbed Acceleration
     acc_mag =  sqrt(acc_x_1.^2+acc_y_1.^2+acc_z_1.^2);
     plot(acc_mag);
     set(gca,'FontSize',15);
     hold on;
      % Magnitude of Perturbed Acceleration
     acc_pert_mag = sqrt(acc_x_2.^2+acc_y_2.^2+acc_z_1.^2);
     plot(acc_pert_mag);
     set(gca,'FontSize',15);
     title('Magnitude of the Acceleration Vector for Satellite-',num2str(s));
     xlabel('Time(seconds)'), ylabel('Acceleration (metres/second^2)');
     legend('Ideal acceleration magnitude','Perturbed acceleration magnitude');     
     xlim([0 10000])
    
    
end


