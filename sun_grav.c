
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
    double pi = 22/7;
    double i_earth = 23.43 ; // degrees
    
    double i_earth_rad = 23.43*pi/180 ; // radians
    
    double rho_sun = 149597870000 ; // meters
    double gamma_sun = 3.9e-14 ; // per s^2
    
    double theta_sat = atan( pos[1]/pos[0]); // radians
    int i,j,k;
    double mul_1[3][3], mul_2[3][1];
    
    // transformation matrix along x axis inclination or earth in solar orbit
    double transformation_matrix_x[3][3] = {
    {1, 0, 0},
    {0, cos(i_earth_rad), sin(i_earth_rad)},
    {0, -sin(i_earth_rad), cos(i_earth_rad)}
    };
    
    double transformation_matrix_x_inverse[3][3]= {
        {1, 0, 0},
        {0, cos(i_earth_rad), -sin(i_earth_rad) },
        {0, sin(i_earth_rad), cos(i_earth_rad) }
    };
    
    // second transformation is around z axis of the satellite
    double transformation_matrix_z[3][3] = {
    { cos(-theta_sat), 0, -sin(-theta_sat)},
    {0,1, 0},
    {sin(-theta_sat), 0, cos(-theta_sat)}
    };
    
    double transformation_matrix_z_inverse[3][3]= {
        { cos(-theta_sat), 0, sin(-theta_sat)},
        {0,1, 0},
        {-sin(-theta_sat), 0, cos(-theta_sat)}
    };
    
        
    // distance between Sun and satellite
    double rho_sun_vector[3][1] = {rho_sun, 0, 0};
    
    // Calculate multiplication of matrices
    
    // transformation_matrix_x_inverse X transformation_matrix_z_inverse = mul_1    
    for(i=0;i< 3;i++)
        {
            for(j=0;j< 3;j++)
            {
                mul_1[i][j] = 0;
                for(k=0;k< 3;k++)
                    {
                        mul_1[i][j] = mul_1[i][j] + transformation_matrix_x_inverse[i][k]*transformation_matrix_z_inverse[j][k];
                    }
            }
            
         }
    
    // mul_1 X rho_sun_vector[3][1] = mul_2[3][1]    
    for ( i = 0; i < 3; ++i) 
    {
        for ( j = 0; j < 1; ++j)
        {
            mul_2[i][j] = 0; // setting the elements of the multiplied vector to zero
        }
    }
    
    for ( i = 0; i < 3; ++i)
    {
        for ( j = 0; j < 1; ++j) 
        {
            for ( k = 0; k < 3; ++k)
            {
                mul_2[i][j] += mul_1[i][k] * rho_sun_vector[k][j];
            }
        }
    }
    
    // mul_2[3][1] ~ mul_2[x,y,z] gives us the elements for further calculation
    double mul_2_x = mul_2[0][1];
    double mul_2_y = mul_2[1][1];
    double mul_2_z = mul_2[2][1];       

    acc_sun_grav[0]= gamma_sun *  mul_2_x; // m/s^2
    acc_sun_grav[1]= gamma_sun *  mul_2_y; // m/s^2
	acc_sun_grav[2]= gamma_sun *  mul_2_z; // m/s^2 
         
}
