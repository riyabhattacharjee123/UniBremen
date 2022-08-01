
/*                                                                    */
/*            Calculate perturbation due to oblateness of Earth
              Adding J2 effects*/
/*                                                                    */
/* ------------------------------------------------------------------ */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "J2_effect.h"

void J2_effect(double *pos,double *acc_J2)
{
    double R;
    double mu;
    double r_dot_dot[3];
    double r_dot [3];
    double r[3];
    
    double radius_earth; 
    double J2; 
    double mass_sat;  
    double sat_length; 
    double sat_breadth; 
    double sat_height; 
    
    mu = 3.986e14;
    radius_earth = 6371000 ; //metres
    J2 = 0.00108263 ; // constant
    mass_sat = 487 ; // Kilograms 
    sat_length = 1.942 ; // metres
    sat_breadth = 3.123 ; // metres
    sat_height = 0.72 ; // metres
    int i;       
    
    R = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]) ; // meters
    
    // calculate the perturbation acceleration in ECI
    acc_J2[0]=3*J2*mu*radius_earth*radius_earth/(2 * pow(R,7))*(5*pow(pos[2],2)-pow(R,2)*pos[0]);
    acc_J2[1]=3*J2*mu*radius_earth*radius_earth/(2 * pow(R,7))*(5*pow(pos[2],2)-pow(R,2)*pos[1]);
	acc_J2[2]=3*J2*mu*radius_earth*radius_earth/(2 * pow(R,7))*(5*pow(pos[2],2)-3*pow(R,2)*pos[2]);

 
         
}
