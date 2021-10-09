
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
    double pi = 22/7;
    double i_earth = 23.43 ; // degrees
    double i_moon = 5.3 ; // degrees 
    double i_earth_rad = 23.43*pi/180 ; // radians
    double i_moon_rad = 5.3*pi/180 ; // radians
    double rho_moon_sat = 384400000 ; // meters
    double gamma_moon = 8.6e-14 ; // per second^2   
    double theta_moon = 2*pi/(2.419e+6); // rad/s moon's rotation angle integrated for 28 days. radians
    int i,j,k;
    double mul_1[3][3], mul_2[3][3], mul_3[3][1];
    
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
    
    // transformation matrix along trasformation y axis of moon inclination
    double transformation_matrix_y2[3][3] = {
    { cos(i_moon_rad), 0, -sin(i_moon_rad)},
    {0,1, 0},
    {sin(i_moon_rad), 0, cos(i_moon_rad)}
    };
    
    double transformation_matrix_y2_inverse[3][3]= {
        { cos(i_moon_rad), 0, sin(i_moon_rad)},
        {0,1, 0},
        {-sin(i_moon_rad), 0, cos(i_moon_rad)}
    };
    
    // transformation around the z-axis of the mooon's rotation angle
    double transformation_matrix_z3[3][3] = {
        {cos(theta_moon), sin(theta_moon), 0},
        {-sin(theta_moon), cos(theta_moon), 0},
        {0, 0, 1}
    };
    
    double transformation_matrix_z3_inverse[3][3] = {
        {cos(theta_moon), -sin(theta_moon), 0},
        {sin(theta_moon), cos(theta_moon), 0},
        {0, 0, 1}
    };
    
    // distance between moon and satellite
    double rho_moon_sat_vector[3][1] = {rho_moon_sat, 0, 0};
    
    // Calculate multiplication of matrices
    
    // transformation_matrix_x_inverse X transformation_matrix_y2_inverse = mul_1    
    for(i=0;i< 3;i++)
        {
            for(j=0;j< 3;j++)
            {
                mul_1[i][j] = 0;
                for(k=0;k< 3;k++)
                    {
                        mul_1[i][j] = mul_1[i][j] + transformation_matrix_x_inverse[i][k]*transformation_matrix_y2_inverse[j][k];
                    }
            }
            
         }
    
    // mul_1 X transformation_matrix_z3_inverse = mul_2
    for(i=0;i< 3;i++)
        {
            for(j=0;j< 3;j++)
            {
                mul_2[i][j] = 0;
                for(k=0;k< 3;k++)
                    {
                        mul_2[i][j] = mul_2[i][j] + mul_1[i][k]*transformation_matrix_z3_inverse[j][k];
                    }
            }
            
         }
    
    // mul_2[3][3] X rho_moon_sat_vector[3][1] = mul_3[3][1]
    for ( i = 0; i < 3; ++i) 
    {
        for ( j = 0; j < 1; ++j)
        {
            mul_3[i][j] = 0; // setting the elements of the multiplied vector to zero
        }
    }
    
    for ( i = 0; i < 3; ++i)
    {
        for ( j = 0; j < 1; ++j) 
        {
            for ( k = 0; k < 3; ++k)
            {
                mul_3[i][j] += mul_2[i][k] * rho_moon_sat_vector[k][j];
            }
        }
    }
    
    // mul_3[3][1] ~ mul_3[x,y,z] gives us the elements for further calculation
    double mul_3_x = mul_3[0][1];
    double mul_3_y = mul_3[1][1];
    double mul_3_z = mul_3[2][1];
    
    // Calculating the RAAN of Moon and the Argument of perigee of the Satellite
    double phi_moon = acos(( mul_3_z/rho_moon_sat)); // radians
    double theta = atan((mul_3_y/mul_3_x)); // radians   

    acc_moon_grav[0]= gamma_moon * sin(phi_moon) * cos(theta); // m/s^2
    acc_moon_grav[1]= gamma_moon * sin(phi_moon) * sin(theta); // m/s^2
	acc_moon_grav[2]=gamma_moon * cos(phi_moon); // m/s^2 
         
}
