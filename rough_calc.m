% Third Body Perturbations due to Moon

i_earth = 23.43 ; % degrees
i_earth_rad = i_earth*pi/180; % rad
i_moon = 5.3 ; % degrees
i_moon_rad = i_moon*pi/180 ; % rad
theta_moon = 2*pi/(2.419e+6); % rad/s moon's rotation angle integrated for 28 days. radians


% First Transformation around x_axis in ECI

transformation_matrix_x = [ 1 0 0;
                            0 cos(i_earth_rad) sin(i_earth_rad); 
                            0 -sin(i_earth_rad) cos(i_earth_rad)];
fprintf('\n transformation matrix =  \n')
disp(transformation_matrix_x)

transformation_matrix_x_inverse = transformation_matrix_x^-1; %inv(transformation_matrix_x)

fprintf('\n transformation matrix inverse =  \n')
disp(transformation_matrix_x_inverse)

transformation_matrix_y2 = [ cos(i_moon_rad) 0 -sin(i_moon_rad);
                             0 1 0;
                             sin(i_moon_rad) 0  cos(i_moon_rad)];
                         
fprintf('\n transformation matrix y2 =  \n')
disp(transformation_matrix_y2)

transformation_matrix_y2_inverse = transformation_matrix_y2^-1; %inv(transformation_matrix_x)

fprintf('\n transformation matrix y2 inverse =  \n')
disp(transformation_matrix_y2_inverse)     

transformation_matrix_z3 = [ cos(theta_moon) sin(theta_moon) 0;
                            -sin(theta_moon) cos(theta_moon) 0;
                            0 0 1];
fprintf('\n transformation matrix z3 =  \n')
disp(transformation_matrix_z3)

transformation_matrix_z3_inverse = transformation_matrix_z3^-1; %inv(transformation_matrix_x)

fprintf('\n transformation matrix z3 inverse =  \n')
disp(transformation_matrix_z3_inverse) 



x_position = 1.1374e+07 ; % meters
y_position = 1.0000e+04 ; % meters
z_position = 8.6603e+03 ; % meters
theta_sat = atan( y_position/x_position);

r_vector_2_sun = [ cos(-theta_sat) 0 -sin(-theta_sat);
                   0 1 0;
                   sin(-theta_sat) 0 cos(-theta_sat)];

r_vector_2_sun_inverse = inv(r_vector_2_sun); 
disp(r_vector_2_sun_inverse);
                        
