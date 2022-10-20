
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
    double lambda = mean_longitude + M + (6892/3600)*sin(M*pi/180) + (72/3600)*sin(2*M*pi/180); //degree    
    // Sun's ecliptic distance
    // distance between the Sun and Earth       
    double rho_sun = 149.619 - 2.499*cos(M*pi/180) - 0.021 * cos(2*M*pi/180) * 1000000000 ; // meters    
    double Mass_sun = 1.989e+33 ; // grams
    // Heliocentric Gravitational constant of Sun = G_sun*Mass_sun
    double mu_sun = 1.32712440042e+20 ; // m^3/s^2
    
    // distance of satellite and Sun in ECI = mean distance of sun and Earth
    // gamma_sun = mu_sun/(geocentric distance of Sun )
    double gamma_sun = 3.9e-14 ; // per s^2   
    // calculate the perturbation acceleration in ECI
    acc_sun_grav[0]= -gamma_sun*rho_sun * cos(lambda*pi/180); // m/s^2
    acc_sun_grav[1]= -gamma_sun*rho_sun * sin(lambda*pi/180) * cos(epsilon*pi/180); // m/s^2
	acc_sun_grav[2]= -gamma_sun*rho_sun * sin(lambda*pi/180) * sin(epsilon*pi/180); // m/s^2 
         
}
