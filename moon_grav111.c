
/*                                                                    */
/*    Calculate perturbation acceleration due to Gravitation of Moon  */
/*                                                                    */
/* ------------------------------------------------------------------ */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "moon_grav.h"

void moon_grav(double *pos,double *acc_moon_grav)
{   
    double pi = 22.0/7;
    double val = pi/180;
    // Julian date fixed for our calculations 01.November.2021 00H:00M:00S
    double JD = 2459519.50000;
    double T = (JD-2451545.0)/36525.0 ;    
    // Mean longitude of moon
    double L_0 = 218.31617 + 481267.88088*T - 1.3972*T ; // degrees    
    //Mean anomaly of moon
    double l = 134.96292 + 477198.86753 * T ; // degrees
    // Mean anomaly of sun
    double l_sun = 357.52543 + 35999.04944 * T ; // degrees
    // Mean angular distance between the Mean and the ascending node 
    double F = 98.272823 + 483202.01873 * T ; // degrees
    // difference between the mean longitudes of Sun and Moon
    double D = 297.85027 + 445267.11135 * T ; // degrees
    
    // Moon's ecliptic longitude
    double lambda = L_0 + 22640/3600 * sin(l*val)+ 769/3600 * sin(2*l*val) \
            -4586/3600 * sin((l-2*D)*val) + 2370/3600 * sin(2*D*val) \
            -668/3600 * sin(l_sun*val) - 412/3600 * sin(2*F*val) \
            -212/3600 * sin(2*(l-D)*val) - 206/3600 * sin((l+l_sun-2*D)*val)\
            +192/3600 * sin((l+2*D)*val) - 165/3600 * sin((l_sun-2*D)*val) \
            + 148/3600 * sin((l-l_sun)*val) - 125/3600 * sin(D*val)\
            - 110/3600 * sin((l+l_sun)*val) - 55/3600 * sin((2*F-2*D)*val) ;
    
    // Moon's latitude
    double beta = 18520/3600 * sin(val*(F+lambda-L_0+ 412/3600 * sin(2*F*val) + 541/3600 * sin(l_sun*val)))\
            -526/3600 * sin((F-2*D)*val) + 44/3600 * sin((l+F-2*D)*val)\
            -31/3600 * sin((-l+F-2*D)*val) - 25/3600 * sin((-2*l+F)*val)\
            -23/3600 * sin((l_sun+F-2*D)*val) + 21/3600 * sin((-l+F)*val)\
            + 11/3600 * sin((-l_sun+F-2*D)*val);
    // Moon's obliquity of the ecliptic = epsilon
    double epsilon = (23+27/60+8/3600) + (46.845/3600 * T)\
            - (0.0059/3600 * pow(T,2)) + ( 0.00181/3600 * pow(T,3)); // degrees
    // Moon's distance from the center of the Earth
    double rho_m = 385000 - 20905*cos(l*val) - 3699* cos((2*D-l)*val)\
            -2956 * cos(2*D*val) - 570 * cos(2*l*val) + 246 * cos(2*(l-D)*val)\
            -205*cos((l_sun-2*D)*val) - 171 * cos((l+2*D)*val) - 152*cos((l+l_sun-2*D)*val)  ; // Km
    double rho_moon = 1000*rho_m ; // m
    // distance of satellite and Moon in ECI = mean distance of Moon and Earth
    // gamma_moon = mu_moon/(geocentric distance of moon )     
    double gamma_moon = 8.6e-14 ; // per second^2       
    // calculate the perturbation acceleration in ECI
    acc_moon_grav[0]= -gamma_moon * rho_moon* cos(lambda*val) * cos(beta*val); // m/s^2
    acc_moon_grav[1]= -gamma_moon * rho_moon* sin(lambda*val) * cos(beta*val) * cos(epsilon*val); // m/s^2
	acc_moon_grav[2]= -gamma_moon * rho_moon* sin(beta*val) * sin(epsilon*val); // m/s^2          
}
