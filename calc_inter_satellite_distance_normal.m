Nr = 64;
k = 1:100 ;
r1= r_start_sat_gs_spherical(1,1); %;600000 %8.897389705668777e+06 % 8.897714448631750e+06; %9.554100016146174e+06;%meters 
r2= r_start_sat_gs_spherical(1,4);%600200

theta_inc_l=[];
theta_inc_j=[];
theta_inc_i = r_start_sat_gs_spherical(1,2) ; %pi/6;


theta_el_l=[];
theta_el_j=[];
theta_el_i = pi/2 - theta_inc_i; %0.963645591331114 ;% radians 0.667321259360872;

phi_i = r_start_sat_gs_spherical(1,3) ; %) pi/3; %-0.681495497830870 %0.542473753252491; %-0.963735944517760 ; % radians
phi_l=[];
phi_j=[];

for ck=1:length(k)    
    if not(mod(ck,Nrx)== 0)
        theta_el_l(1,ck) = acos( cos(theta_el_i) + (2*ck/Nrx) ); % radians
        theta_el_j(1,ck) = acos( cos(theta_el_l(1,ck)) + (2*ck/Nrx) ); % radians
        
        theta_inc_l(1,ck) = pi/2 - theta_el_l(1,ck) ; %radians
        theta_inc_j(1,ck) = pi/2 - theta_el_j(1,ck) , % radians
        
        %phi_l(1,ck) = acos(((2*ck/Nrx) + ...
         %                    cos(phi_i)* sin(theta_el_i))...
          %                 /cos(theta_el_l(1,ck)) ); % radians
        phi_l(1,ck) = phi_i; % acos(cos(phi_i)+(2*ck/Nrx));
        phi_j(1,ck) = phi_i;% acos(cos(phi_l(1,ck))+(2*ck/Nrx));
        %phi_j(1,ck) = acos(((2*ck/Nrx) + ...
         %                   cos(phi_l(1,ck))* sin(theta_el_l(1,ck)))...
          %                /cos(theta_el_j(1,ck)) ); % radians                         
        
    end
end

%last = 0.472150471443591 - 0.881517970500656i ;
 %8.897347307200497e+06 %8.896866471827822e+06; %8.285765311388365e+06;%meters 
r3= 8.897380839692418e+06 %8.897537133897793e+06; %9.303186844548663e+06;%meters 

distance_1= [];
distance_2= [];
for cdist=1:length(k)
   %distance_2(1,cdist) = sqrt(r1^2 + r2^2 - 2*r1*r2...
    %                          *(cosd((theta_i-theta_l(1,cdist))*180/pi))); %meters
    
    
                          
   distance_2(1,cdist) = sqrt(r1^2 + r2^2 - 2*r1*r2...
                              *( sind(theta_el_i*180/pi)...
                              * sind(theta_el_l(1,cdist)*180/pi)...
                              *( cosd((phi_i-phi_l(1,cdist))*180/pi))...
                              + cosd(theta_el_i*180/pi)...
                              *cosd((theta_el_l(1,cdist))*180/pi) )); %meters
                          
  x_1 = r1*sind(180/pi * (theta_inc_i))*sind(180/pi * ((2*pi)-phi_i)) ;
  x_2 = r1*sind(180/pi * (theta_inc_l(1,cdist)))*sind(180/pi * ((2*pi)-phi_l(1,cdist))) ;
  
  y_1 = r1*cosd(180/pi * ((pi/2)-phi_i));
  y_2 = r2*cosd(180/pi * ((pi/2)-phi_l(1,cdist)));

  z_1 = r1*cosd(180/pi * theta_inc_i)*sind(180/pi * ((2*pi)-phi_i)) ;
  z_2 = r2*cosd(180/pi * theta_inc_l(1,cdist))*sind(180/pi * ((2*pi)-phi_l(1,cdist))) ;

  distance_1(1,cdist) = sqrt((x_1-x_2)^2+(y_1-y_2)^2+(z_1-z_2)^2);
                          
end

