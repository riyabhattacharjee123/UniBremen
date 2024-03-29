% ref : https://core.ac.uk/download/pdf/4268217.pdf

% consider a Leader-Chaser formation. ONE Chaser. One Leader exist in the
% same orbital plane.
% The Chaser move in circular orbit around the leader. The leader is
% also in a circular orbit around the Earth.

% The below code calculates the initial position and velocities of Chaser
% --in Clohessy-Wiltshire Equations in General Circular Orbit

%clear all
% rho = Relative distance between the Chaser and the leader in ECI-y axis
rho = side ; % metres  %20000; %44000; 
mu = 3.9857 *10^(14); % constant
disp('rho (m)');
disp(rho) ;
radius_earth = 6371000 ; % metres

% Initializing the starting position of leader satellite in ECI
x_leader_altitude = 370000; %metres
x_leader = radius_earth + x_leader_altitude; %metres
y_leader = 0; %metres
z_leader = 0 ; %metres

% calculate the initial velocity of the leader satellite in circular orbit
% using vis-viva equation
% since the motion is only in y-direction n ECI. we consider only x-y plane
xdot_leader = 0; % metres/second
ydot_leader = sqrt ( mu / x_leader); % metres/second
zdot_leader = 0; % metres/second

% Calculate the initial state vectr of the chaser satellites using the CW
% --equations in the LVLH frame.
% We initialize the time, orbital angles, integration constants
% --mean motion
t = 0 ; %seconds
alpha = deg2rad(120 + [0 120].'); % radians, 180 degrees for the chaser
beta = deg2rad(15) ; % radians , 15 degrees

% Defining the integration constants of the HCW equations
c1 = rho ;
c2 = (sqrt(3)/2)*rho;
c3 = 0;

n = sqrt(mu/(x_leader^3)); % mean motion of target satellite 
theta_dot = [0 0 n]; % rad/second angular velocity of LVLH frame wrt ECI around z-axis

% Calculate the position of chasers in LVLH frame using the CW formulae
% Add the translation term to get the values in ECI
x_0 = (c1/2) * sin(n*t + alpha); % metres in LVLH-x
x_c = x_leader + x_0; %metres in ECI-x

y_0 = (c1) * cos( n*t + alpha ) + c3 ; % metres in LVLH-y
y_c = y_0 + y_leader ; %metres in ECI-y

z_0 = [c2 * sin(n*t + beta ) ; c2 * sin(n*t + beta )] ; % metres in LVLH-z
z_c = z_0 + z_leader ; %metres in ECI-z

% Display the above calculated values
disp('absolute_radius_chaser x_C');
disp(x_0);
disp(x_c);
disp('absolute_radius_chaser y_C');
disp(y_0);
disp(y_c);
disp('absolute_radius_chaser z_C');
disp(z_0);
disp(z_c);

% Calculate the initial velocities of the chasers in LVLH frame using
% --formulae in CW
% Translate the result into ECI frame

rho_eci = [x_0 y_0 z_0]; % metres distance from leader to chaser in LVLH. Matrix format

disp ('rho_eci')
disp (rho_eci)
disp('repmat theta_dot')
disp(repmat(theta_dot,1,1))

rotation_term = cross(repmat(theta_dot,2,1),rho_eci,2); % we get the rotation translation matrix
xdot_rotation = rotation_term(:,1) ;
ydot_rotation = rotation_term(:,2) ;
zdot_rotation = rotation_term(:,3) ;


xdot_0 = (c1/2) * n * cos (alpha); % metres/second in LVLH-x
xdot_c = xdot_leader + xdot_0 + xdot_rotation ; % metres/second in ECI-x
disp('xdot_c');
disp(xdot_c);

ydot_0 = (-c1) * n * sin(alpha); % metres/second in LVLH-y
ydot_c = ydot_leader + ydot_0 + ydot_rotation ; % metres/second in ECI-y
disp('ydot_c');
disp(ydot_c);

zdot_0 = c2 * n * cos(beta); % metres/second in LVLH-z
zdot_c = zdot_leader + zdot_0 + zdot_rotation ; % metres/second in ECI-z
disp('zdot_c');
disp(zdot_c);

