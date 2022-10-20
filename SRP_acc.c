#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "SRP_acc.h"

void SRP_acc(double *pos,double *SRP_acc_out)
{
    long double pi = 22.0/7.0;      
    double R_sun_earth = 149633425000; //meters
    //assuming same because of very small difference
    double R_sat_sun = R_sun_earth;    
    double S_0=1368; // Watt/m^2
    double mass_sat = 487000 ; // grams 
    double sat_length = 1.942 ; // metres
    double sat_breadth = 3.123 ; // metres
    double sat_height = 0.72 ; // metres
    double Area_sat = sat_length * sat_breadth ; // meters^2
    double c=3e8; // m/s
    double alpha=0.4; // reflection coeff
    double gamma=0;
    // Julian date fixed for our calculations 01.November.2021 00H:00M:00S
    double JD=2459519.50000;   
    double n=JD-2451545;
    double L= 280.460+0.9856474*n; // degrees
    double g=357.528+0.9856003*n; //degrees
    double theta=L+1.915*sin(g*pi/180)+0.020*sin(2*g*pi/180);//degrees
    double epsilon=23.439-0.0000004*n;//degrees
    // calculate the perturbation acceleration in ECI
    SRP_acc_out[0]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*\
            S_0/c*(1+alpha)*pow(cos(gamma),2))*cos(theta*pi/180);//m/s^2
     SRP_acc_out[1]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*\
             S_0/c*(1+alpha)*pow(cos(gamma),2))*cos(epsilon*pi/180)\
             *sin(theta*pi/180);//m/s^2
      SRP_acc_out[2]=-pow((R_sun_earth/R_sat_sun),2)*(Area_sat/mass_sat*\
              S_0/c*(1+alpha)*pow(cos(gamma),2))*sin(epsilon*pi/180)\
              *sin(theta*pi/180);//m/s^2
}
