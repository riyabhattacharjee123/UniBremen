
/*                                                                    */
/* Calculate perturbation acceleration due to Gravitation of the Sun  */
/*                                                                    */
/* ------------------------------------------------------------------ */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "sun_grav.h"

void sun_grav(double *pos,double *acc_sun_grav)
{   
    double pi = 22./7;  
    double val = pi/180;
    double rho_sun_sat[3];
    double sun_cood[3], part_one[3], part_two[3];
    double mu_earth = 3.9857e+14;
    double mu_sun = 1.32712440042e+20 ; // m^3/s^2
    // Julian date fixed for our calculations 01.November.2021 00H:00M:00S
    double JD = 2459519.50000;
    double T = (JD-2451545.0)/36525.0 ;
    //obliquity of the ecliptic. inclination of ecliptic relative to earth's equator
    double epsilon = 23.43 ; // degrees      
    //Mean longitude of sun (omega + OMEGA)
    double mean_longitude = 282.94000 ;// degrees
    // Mean anomaly of Sun. Pretending Sun orbits around Earth
    double M = 357.5256 + 35999.049 * T ; // degree
    
    // Sun's ecliptic longitude = lambda
    double lambda = mean_longitude + M + (6892.0/3600.0)*sin(M*val)\
            + (72.0/3600.0)*sin(2*M*val); //degree    
    // Sun's ecliptic distance
    // distance between the Sun and Earth       
    double rho_sun = (149.619 - 2.499*cos(M*val)\
            - 0.021 * cos(2*M*val)) * 1000000000 ; // meters    
    
    sun_cood[0] = rho_sun* cos(lambda*val);
    sun_cood[1] = rho_sun* sin(lambda*val) * cos(epsilon*val);
    sun_cood[2] = rho_sun* sin(lambda*val) * sin(epsilon*val);    
    
    // Heliocentric Gravitational constant of Sun = G_sun*Mass_sun
    
    rho_sun_sat[0] = sun_cood[0] - pos[0];
    rho_sun_sat[1] = sun_cood[1] - pos[1];
    rho_sun_sat[2] = sun_cood[2] - pos[2];
    double dr_a = pow((pow(rho_sun_sat[0],2) + pow(rho_sun_sat[1],2)+ pow(rho_sun_sat[2],2)),1.5);
    //double dr_a = pow(mag_rho_sun_sat,3);
    double dr_b = pow((pow(sun_cood[0],2) + pow(sun_cood[1],2)+ pow(sun_cood[2],2)),1.5);
    //double dr_b = pow(mag_rho_sun,3);
       
    double mag_sat_pos = sqrt(pow(pos[0],2) + pow(pos[1],2)+ pow(pos[2],2));    
    
    part_one[0] = (mu_sun/dr_a)* rho_sun_sat[0];
     part_one[1] = (mu_sun/dr_a)* rho_sun_sat[1];
      part_one[2] = (mu_sun/dr_a)* rho_sun_sat[2];
      
    part_two[0] = (mu_sun/dr_b)*sun_cood[0];
     part_two[1] = (mu_sun/dr_b)*sun_cood[1];
      part_two[2] = (mu_sun/dr_b)*sun_cood[2];
    
    
    // calculate the perturbation acceleration in ECI
    //acc_sun_grav[0]=  mu_sun *(  (rho_sun_sat[0]/dr_a) - (sun_cood[0]/dr_b) ); // m/s^2
    //acc_sun_grav[1]=  mu_sun * ((rho_sun_sat[1]/dr_a)- (sun_cood[1]/dr_b)); // m/s^2
	//acc_sun_grav[2]=  mu_sun *((rho_sun_sat[2]/dr_a)- (sun_cood[2]/dr_b) ); // m/s^2 
     acc_sun_grav[0]=   part_one[0]-part_two[0];
    acc_sun_grav[1]=  part_one[1]-part_two[1];
    acc_sun_grav[2]=  part_one[2]-part_two[2];
               
}
