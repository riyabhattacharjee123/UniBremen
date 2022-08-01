
/*                                                                    */
/*            Calculate perturbation due to air drag force
/*                                                                    */
/* ------------------------------------------------------------------ */ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "air_drag.h"

void air_drag(double *pos,double *vel,double *acc_drag)
{
    double R;
    double mu;
    double h;
    double C_D;
    double F10;
    double Ap;
    double radius_earth;
    double mass_sat;  
    double sat_length; 
    double sat_breadth; 
    double sat_height; 
    double mu_air;
    double T;
    double H;
    double rho_air;
    double Area_sat;
    double air_drag_acc;
    double cross_product[3];
    double vel_sat_rel[3];
    radius_earth = 6371000 ; //metres
    C_D = 2.2 ; // Drag coefficient
    h = 500 ; // Kilometres
    // data given from http://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_nowcast.txt
    // date 01 November 2021
    F10 = 97.7 ;
    Ap = 8;
    
    mu = 3.986e14;
   
    mass_sat = 487 ; // Kilograms 
    sat_length = 1.942 ; // metres
    sat_breadth = 3.123 ; // metres
    sat_height = 0.72 ; // metres
    int i;
    
    T = 900 + 2.5 *(F10-70) + 1.5* Ap;
    mu_air = 27 - 0.012 * (h-200);
    H = T / mu_air ;
    
    double omega_earth_avg =  7.2921159e-5 ; // rad/second
    
    //https://www.spaceacademy.net.au/watch/debris/atmosmod.htm
    //https://www.researchgate.net/publication/278327611_The_Effect_of_Aerodynamic_Drag_Forces_on_the_Formation_Flying_of_Satellites
    rho_air = 6e-10 * exp (-(h-175)/H); // Kg/m^3
    Area_sat = sat_length * sat_breadth ;
    
    // Calculate cross product of omega_earth and inertial position of satellite

    cross_product[0] = -omega_earth_avg*pos[1];
    cross_product[1] = -omega_earth_avg*pos[0];
    cross_product[2] = 0;
    
    // velocity_spacecraft_relative = velocity_spacecraft - omega_earth X position_spacecraft
    vel_sat_rel[0] = vel[0]-cross_product[0];
    vel_sat_rel[1] = vel[1]-cross_product[1];
    vel_sat_rel[2] = vel[2]-cross_product[2];
    
    // calculate the magnitude of velocity of satellite in inertial frame    
    double vel_rel_mag = sqrt(pow(vel_sat_rel[0],2) + pow(vel_sat_rel[1],2)+ pow(vel_sat_rel[2],2)); // m/s
    // calculate the perturbation acceleration in ECI
    acc_drag[0]= -(1./2) * C_D * rho_air * (Area_sat/mass_sat)* vel_rel_mag * vel_sat_rel[0]  ; // m/second^2 
    acc_drag[1]= -(1./2) * C_D * rho_air * (Area_sat/mass_sat)* vel_rel_mag * vel_sat_rel[1]  ; // m/second^2 
    acc_drag[2]= -(1./2) * C_D * rho_air * (Area_sat/mass_sat)* vel_rel_mag * vel[2]  ; // m/second^2
    
}
