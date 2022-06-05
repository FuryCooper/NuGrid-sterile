#ifndef ZZCNEUTRINO_H
#define ZZCNEUTRINO_H


#include <stdio.h>
#include <gsl/gsl_rng.h>

/*double xi_expan = 0.5; 						//neutrino chemical potential over temperature
double T_neu0 = 1.9;					//neutrino temperature at z=0 (Kalvein)
double mass_neu_expan = 0.1;				//neutrino mass (ev)

double Gr = 6.67384E-11;   	//m^3*kg^-1*s^-2
double kb = 1.3806488E-23;	//m^2*kg*s^-2*k^-1
double c_light = 299792458.0; 		//m*s^-1
double cc = 299792458.0*299792458.0;
double ktoev = 1.0/11605.0;
double hbar = 1.05457173E-34;	//m^2*kg*s^-1
double H0 = 2.17489E-18;		//s^-1

*/
extern double xi_expan; 						//neutrino chemical potential over temperature
extern double T_neu0;					//neutrino tempearture at z=0
extern double mass_neu_expan;				//neutrino mass

extern double hbar;
extern double kb;
extern double Gr;
extern double c_light;
extern double cc;
extern double ktoev;
extern double hbar;
extern double H0;
extern double h;
extern double roneu0;
double neutrino_partition(double p, void *param);
double neutrino_integration(double a);

#endif
