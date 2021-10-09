#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "SRP_acc.h"

#define pi 3.142857

void SRP_acc(double *pos,double *SRP_out)
{
    double R_sun_earth;
    double r_dot_dot[3];
    double r_dot [3];
    double r[3];
    double n;
    double R_sat_sun; 
    double mass_sat;  
    double sat_length; 
    double sat_breadth; 
    double sat_height; 
    double Area_sat;
    double S_0;
    double c;
    double alpha;
    double gamma;
    double JD, L, g, theta, epsilon;
    double acc_srp[3];
    Area_sat = sat_length * sat_breadth ;
    R_sun_earth=149633425000; //meters
    R_sat_sun=R_sun_earth; // assuming it same because of very small difference
    
    S_0=1368; // Watt/m^2
    mass_sat = 487000 ; // grams 
    sat_length = 1.942 ; // metres
    sat_breadth = 3.123 ; // metres
    sat_height = 0.72 ; // metres
    c=3e8; // m/s
    alpha=0.4;
    gamma=0;
    Area_sat = sat_length * sat_breadth ;
    JD=21277; // Julian Date 
    
    // Calculate acceleration SRP
    // for loop for each row
         
    n=JD-2451545;
    L= 280.460+0.9856474*n; // degrees
    g=357.528+0.9856003*n; //degrees
    theta=L+1.915*sin(g*pi/180)+0.020*sin(2*g*pi/180);
    epsilon=23.439-0.0000004*n;
    
    acc_srp[0]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*S_0/c*(1+alpha)*pow(cos(gamma),2))*cos(theta*pi/180);
    acc_srp[1]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*S_0/c*(1+alpha)*pow(cos(gamma),2))*cos(epsilon*pi/180)*sin(theta*pi/180);
    acc_srp[2]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*S_0/c*(1+alpha)*pow(cos(gamma),2))*sin(epsilon*pi/180)*sin(theta*pi/180);
}
