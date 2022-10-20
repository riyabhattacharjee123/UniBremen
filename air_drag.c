/* ------------------------------------------------------------------ */ 
/*      Calculate perturbation due to atmospheric drag force          */
/* ------------------------------------------------------------------ */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "air_drag.h"

void air_drag(double *pos,double *vel,double *acc_drag)
{   
    double cross_product[3];
    double vel_sat_rel[3];
    double radius_earth = 6371000 ; //metres
    double C_D = 2.2 ; // Drag coefficient
    double h = 500 ; // Kilometres. Height of satellite
    // data taken from 
    //http://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_nowcast.txt
    // date 01 November 2021
    double F10 = 97.7 ;
    double Ap = 8;
    
    double mu = 3.986e14; // Earth´s gravitational constant   
    double mass_sat = 487 ; // Kilograms 
    double sat_length = 1.942 ; // metres
    double sat_breadth = 3.123 ; // metres
    double sat_height = 0.72 ; // metres 
    
    double T = 900 + 2.5 *(F10-70) + 1.5* Ap;
    double mu_air = 27 - 0.012 * (h-200);
    double H = T / mu_air ;    
    double omega_earth_avg =  7.2921159e-5 ; // rad/second    
    double rho_air = 6e-10 * exp (-(h-175)/H); // Kg/m^3
    double Area_sat = sat_length * sat_breadth ;
    
    // Cross product of omega_earth and inertial position of satellite
    cross_product[0] = -omega_earth_avg*pos[1];
    cross_product[1] = -omega_earth_avg*pos[0];
    cross_product[2] = 0;
    
    // velocity_spacecraft_relative 
    // = velocity_spacecraft - omega_earth X position_spacecraft
    vel_sat_rel[0] = vel[0]-cross_product[0];
    vel_sat_rel[1] = vel[1]-cross_product[1];
    vel_sat_rel[2] = vel[2]-cross_product[2];
    
    // calculate the magnitude of velocity of satellite in inertial frame    
    double vel_rel_mag = sqrt(pow(vel_sat_rel[0],2)+pow(vel_sat_rel[1],2)\
            +pow(vel_sat_rel[2],2)); // m/s
    // calculate the perturbation acceleration in ECI
    acc_drag[0]= -(1./2)*C_D*rho_air*(Area_sat/mass_sat)\
            *vel_rel_mag*vel_sat_rel[0]  ; // m/second^2 
     acc_drag[1]= -(1./2)*C_D*rho_air*(Area_sat/mass_sat)\
             *vel_rel_mag*vel_sat_rel[1]  ; // m/second^2 
      acc_drag[2]= -(1./2)*C_D*rho_air*(Area_sat/mass_sat)\
              *vel_rel_mag*vel_sat_rel[2]  ; // m/second^2    
}
