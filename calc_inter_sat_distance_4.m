Nr = 64;
k = 1:100 ;
theta_m=[];
theta_l=[];
theta_j=[];
theta_i = 0.958120924147366;% radians 0.667321259360872;

phi_i = 0.549625418639461 ; % radians
phi_l=[];
phi_j=[];
phi_m=[];

for ck=1:length(k)    
    if not(mod(ck,Nrx)== 0)
        theta_l(1,ck) = acos( cos(theta_i) + (2*ck/Nrx) ); % radians
        theta_j(1,ck) = acos( cos(theta_l(1,ck)) + (2*ck/Nrx) ); % radians
        theta_m(1,ck) = acos( cos(theta_j(1,ck)) + (2*ck/Nrx) ); % radians
        
        phi_l(1,ck) = asin(((2*ck/Nrx) + ...
                              sin(phi_i)* cos(theta_i))...
                             /cos(theta_l(1,ck)) ); % radians
        %phi_l(1,ck) = ck/(2*Nrx)  -  theta_i ;
        %phi_j(1,ck) = ck/(2*Nrx)  -  theta_l(1,ck) ;
        phi_j(1,ck) = asin(((2*ck/Nrx) + ...
                             sin(phi_l(1,ck))* cos(theta_l(1,ck)))...
                            /cos(theta_j(1,ck)) ); % radians  
        phi_m(1,ck) = asin(((2*ck/Nrx) + ...
                             sin(phi_j(1,ck))* cos(theta_j(1,ck)))...
                            /cos(theta_m(1,ck)) ); % radians  
        
    end
end

last = 0.472150471443591 - 0.881517970500656i ;
r1= 8.827350075604016e+06 %meters 
r2= 8.827306732242115e+06 %8.896866471827822e+06; %8.285765311388365e+06;%meters 
r3= 8.827293684304256e+06 %8.897537133897793e+06; %9.303186844548663e+06;%meters 
r4= 8.827337027730227e+06

distance_1= [];
%distance_2= [];
for cdist=1:length(k)
   %distance_2(1,cdist) = sqrt(r1^2 + r2^2 - 2*r1*r2...
    %                          *(cosd((theta_i-theta_l(1,cdist))*180/pi))); %meters 
                          
   distance_1(1,cdist) = sqrt(r1^2 + r2^2 - 2*r1*r2...
                              *( sind(theta_i*180/pi)...
                              * sind(theta_l(1,cdist)*180/pi)...
                              *( cosd((phi_i-phi_l(1,cdist))*180/pi))...
                              + cosd(theta_i*180/pi)...
                              *cosd((theta_l(1,cdist))*180/pi) )); %meters
end

