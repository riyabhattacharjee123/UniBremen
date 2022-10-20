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
r_vector = [x_position; y_position; z_position];

r_position_magnitude = sqrt(x_position^2 + y_position^2 + z_position^2) ; % meters

disp('r_position in meters')
disp(r_position_magnitude)


% Step 4
% Third Body Perturbations due to Moon

i_earth = 23.43 ; % degrees
i_moon = 5.3 ; % degrees

% First Transformation around x_axis in ECI

transformation_matrix_x = [ 1 0 0;
                            0 cosd(i_earth) sind(i_earth); 
                            0 -sind(i_earth) cosd(i_earth)];

r_vector_1 = transformation_matrix_x * r_vector;

disp(r_vector_1);

% The second transformation is a rotation of 5.3 degrees around the y2 axis

transformation_matrix_y2 = [ cosd(i_moon) 0 -sind(i_moon);
                             0 1 0;
                             sind(i_moon) 0  cosd(i_moon)];

disp(transformation_matrix_y2);

r_vector_2 = transformation_matrix_y2 * r_vector_1;

disp(r_vector_2);

% Third transformation around the z-axis of the mooon's rotation angle

theta_moon = 2*pi/(2.419e+6); % rad/s moon's rotation angle integrated for 28 days. radians

transformation_matrix_z3 = [ cos(theta_moon) sin(theta_moon) 0;
                            -sin(theta_moon) cos(theta_moon) 0;
                            0 0 1];

rho_moon_sat = 384400000 ;
r_vector_3 = (transformation_matrix_z3)^(-1)* rho_moon_sat ;

disp(r_vector_3);

rrr = inv(transformation_matrix_x) * inv(transformation_matrix_y2) * inv(transformation_matrix_z3) * [rho_moon_sat; 0;0];

disp(rrr);

rrr_x = rrr(1,:);
rrr_y = rrr(2,:);
rrr_z = rrr(3,:);

phi_moon = acos(( rrr_z/rho_moon_sat)); % radians
theta = atan((rrr_y/rrr_x)); % radians

%moon_acc = (1/mass_sat*1000)*[ 

gamma_moon = 8.6e-14 ; % /second sq

acc_moon_x = gamma_moon * sin(phi_moon) * cos(theta);
acc_moon_y = gamma_moon * sin(phi_moon) * sin(theta);
acc_moon_z = gamma_moon * cos(phi_moon);
fprintf('\n Acceleration_moon_gravitation (m/s^2) = [%g %g %g] m/s^2 \n',acc_moon_x, acc_moon_y, acc_moon_z)


% Perturbation accelerations due to Sun

gamma_sun = 3.9e-14 ; % / second square
theta_sat = atan( y_position/x_position);
rho_sun = 149597870000 ; % meters

% first transformation is around x axis

% second transformation is around z axis 

r_vector_2_sun = [ cos(-theta_sat) 0 -sin(-theta_sat);
                   0 1 0;
                   sin(-theta_sat) 0 cos(-theta_sat)]

force_vector_sun = transformation_matrix_x^(-1) * r_vector_2_sun^(-1)*[rho_sun;0;0];

acc_sun_x = (gamma_sun)* force_vector_sun(1,:);
acc_sun_y = (gamma_sun)* force_vector_sun(2,:);
acc_sun_z = (gamma_sun)* force_vector_sun(3,:);

fprintf('\n Acceleration_SUN_gravitation (m/s^2) = [%g %g %g] m/s^2 \n',acc_sun_x, acc_sun_y, acc_sun_z)





