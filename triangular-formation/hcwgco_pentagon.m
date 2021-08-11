% ref : https://core.ac.uk/download/pdf/4268217.pdf

% consider a Pentagon formation. Five Chasers. One Leader that is
% imaginary. and is located at the centre of the equilateral Pentagon
% formed by the Five satellites. This centroid coincides with the origin
% of the LVLH frame.
% The satellites move in circular orbits around the leader. The leader is
% also in a circular orbit around the Earth.

% The below code calculates the initial position and velocities of chasers
% --in Clohessy-Wiltshire Equations in General Circular Orbit

side = 1200 ; % metres , side of the equilateral Pentagon

% rho = Relative distance between the centroid of the Pentagon and the 
% --respective deputies in metres
rho = side / (2 * cosd(54)); % metres

mu = 3.9857 *10^(14); % constant
disp('rho (m)');
disp(rho) ;
radius_earth = 6371000 ; %metres

% Initializing the starting position of leader satellite in ECI
x_leader_altitude = 500000; %metres
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
alpha = 5.18363; % radians, 9 +72  deg for the five phases each time
beta = 1.29154 ; % radians , 74 degrees

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

z_0 = c2 * sin(n*t + alpha ) ; % metres in LVLH-z
z_c = z_0 + z_leader ; %metres in ECI-z

% Display the above calculated values
disp('absolute_radius_chaser x_C');
disp(x_c);
disp('absolute_radius_chaser y_C');
disp(y_c);
disp('absolute_radius_chaser z_C');
disp(z_c);

% Calculate the initial velocities of the chasers in LVLH frame using
% --formulae in CW
% Translate the result into ECI frame

rho_eci = [x_0 y_0 z_0]; % metres distance from leader to chaser in LVLH. Matrix format

rotation_term = cross(theta_dot, rho_eci); % we get the rotation translation matrix
xdot_rotation = rotation_term(1,1) ;
ydot_rotation = rotation_term(1,2) ;
zdot_rotation = rotation_term(1,3) ;


xdot_0 = (c1/2) * n * cos (alpha); % metres/second in LVLH-x
xdot_c = xdot_leader + xdot_0 + xdot_rotation ; % metres/second in ECI-x
disp('xdot_c');
disp(xdot_c);

ydot_0 = (-c1) * n * sin(alpha); % metres/second in LVLH-y
ydot_c = ydot_leader + ydot_0 + ydot_rotation ; % metres/second in ECI-y
disp('ydot_c');
disp(ydot_c);

zdot_0 = c2 * n * cos(alpha); % metres/second in LVLH-z
zdot_c = zdot_leader + zdot_0 + zdot_rotation ; % metres/second in ECI-z
disp('zdot_c');
disp(zdot_c);


