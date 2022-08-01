% step 1
% Calculate perturbation due to oblateness of Earth
% Adding J2 effects

mu = 3.9857 *10^(14); % constant
radius_earth = 6371000 ; %metres
J2 = 0.00108263 ; % constant
mass_sat = 487 ; % Kilograms 
sat_length = 1.942 ; % metres
sat_breadth = 3.123 ; % metres
sat_height = 0.72 ; % metres

x_position = 1.1374e+07 ; % meters
y_position = 1.0000e+04 ; % meters
z_position = 8.6603e+03 ; % meters

r_position_magnitude = sqrt(x_position^2 + y_position^2 + z_position^2) ; % meters

disp('r_position in meters')
disp(r_position_magnitude)

part_1 = -mu/(r_position_magnitude)^3 ;
part_2 = (3/2) * J2 * (radius_earth/r_position_magnitude)^2 ;
part_3 = 5 * (z_position/r_position_magnitude)^2 -1 ;
part_4 = 1 - part_2 * part_3 ;
part_5 = part_1 * part_4 ;

x_acc_J2 = x_position * part_5 ;
y_acc_J2 = y_position * part_5 ;
z_acc_J2 = z_position * part_5 ;


% step 2
% Acceleration due to Atmospheric Drag

C_D = 2.2 ; % Drag coefficient

h = 500 ; % Kilometres
F10 = (70+300)/2 ;
Ap = 4 ; % http://www-app3.gfz-potsdam.de/kp_index/qlyymm.html
T = 900 + 2.5 *(F10-70) + 1.5* Ap;
mu_air = 27 - 0.012 * (h-200);
H = T / mu_air ;
%https://www.spaceacademy.net.au/watch/debris/atmosmod.htm
rho_air = 6e-10 * exp (-(h-175)/H); % Kg/m^3

velocity_sat = 5.9189e+03 ; % metres/second
Area_sat = sat_length * sat_breadth ;

air_drag_acc_y = -(1/2) * C_D * rho_air * velocity_sat^2 * (Area_sat/mass_sat); % m/second^2

disp('accelerations m/s2')
x_acc = x_acc_J2 ;
y_acc = y_acc_J2 + air_drag_acc_y ;
z_acc = z_acc_J2 ;

disp(x_acc);
disp(y_acc);
disp(z_acc);


