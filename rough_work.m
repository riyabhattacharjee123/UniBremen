r1= 8.897389705668777e+06 ;
r2= 8.897347307200495e+06 ;

theta_el_1 = 2.534395591331115; % pi/2 - (-0.963645591331114); % radians
theta_deg_el = 180/pi*theta_el_1 ;
theta_el_2 = acos( cos(theta_el_1) + (2*1/8) ); % radians

theta_inc_1 = -0.963645591331114; %pi/2 - theta_el_1 ; % radians
theta_deg_inc = 180/pi * theta_inc_1;
theta_inc_2 = pi/2 - theta_el_2 ; % radians

phi_1 =  0.542473753252491 ; % radians
phi_2 = acos(cos(phi_1)+(2*1/8)); % radians

% distance between two points

x_1 = r1*sind(180/pi * (theta_inc_1))*sind(180/pi * ((2*pi)-phi_1)) ;
x_2 = r2*sind(180/pi * theta_inc_2)*sind(180/pi * ((2*pi)-phi_2)) ;

y_1 = r1*cosd(180/pi * ((pi/2)-phi_1));
y_2 = r2*cosd(180/pi * ((pi/2)-phi_2));

z_1 = r1*cosd(180/pi * theta_inc_1)*sind(180/pi * ((2*pi)-phi_1)) ;
z_2 = r2*cosd(180/pi * theta_inc_2)*sind(180/pi * ((2*pi)-phi_2)) ;

distance = sqrt((x_1-x_2)^2+(y_1-y_2)^2+(z_1-z_2)^2);



