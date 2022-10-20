/* --------------------------------------------------------------------- */ 
/*  Calculate perturbation due to oblateness of Earth Adding J2 effects  */
/* --------------------------------------------------------------------- */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "J2_effect.h"

void J2_effect(double *pos,double *acc_J2)
{
    double mu = 3.986e14;
    double radius_earth = 6371000 ; //metres
    double J2 = 0.00108263 ; // constant   
    double R = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]); //meters
    
    // calculate the perturbation acceleration in ECI
    acc_J2[0]=( 3*J2*mu*radius_earth*radius_earth ) / (2 * pow(R,7))\
            *( ( 5*pow(pos[2],2)-pow(R,2) )*pos[0] ); // m/s^2
    acc_J2[1]=( 3*J2*mu*radius_earth*radius_earth ) / (2 * pow(R,7))\
            *( ( 5*pow(pos[2],2)-pow(R,2) )*pos[1] ); // m/s^2
	acc_J2[2]=( 3*J2*mu*radius_earth*radius_earth ) / (2 * pow(R,7))\
            *( ( 5*pow(pos[2],2)-3*pow(R,2) )*pos[2] ); // m/s^2         
}
