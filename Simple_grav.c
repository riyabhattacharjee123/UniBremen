/* ------------------------------------------------------------------ */
/*                                                                    */
/*       Simple gravity: Just calculates acceleration by newton gravitational law*/
/*                                                                    */
/* ------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "simple_grav.h"


void simple_grav(double *pos,double *output)
{     
    double r_dot_dot[3];    
    double mu = 3.986e14;   
       
    // Calculate acceleration    
    double radius_square = (pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    r_dot_dot[0] = -mu/radius_square*pos[0]/sqrt(radius_square); 
    r_dot_dot[1] = -mu/radius_square*pos[1]/sqrt(radius_square); 
    r_dot_dot[2] = -mu/radius_square*pos[2]/sqrt(radius_square); 
    
    output[0] = r_dot_dot[0]; 
    output[1] = r_dot_dot[1]; 
    output[2] = r_dot_dot[2];      
}
