/* These amendments zzcneutrino.c and zzcneutrino.h add the influence of massive and degenerate of cosmic
neutrinos into the N-body simulation. A simple test of Hubble constant versus expansion coefficient
relation shows the influence may not be neglectable.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"
#include "zzcneutrino.h"

//static double logTimeBegin;
//static double logTimeMax;
//double T_neu;
//double xi_expan, T_neu0, mass_neu_expan, Gr, kb, c_light, cc, ktoev, hbar, H0;
//double xi_expan = 0.0;
//double xi_expan = 0.6528233;            //neutrino chemical potential over temperature
//double T_neu0 = 1.94501196;          //neutrino temperature at z=0 (Kalvein)
//double mass_neu_expan = 0.0;        //neutrino mass (ev)
/*
double Gr = 6.67384E-11;    //m^3*kg^-1*s^-2
double kb = 1.3806488E-23;  //m^2*kg*s^-2*k^-1
double c_light = 299792458.0;     //m*s^-1
double cc = 299792458.0*299792458.0;
double ktoev = 1.0/11605.0;
double hbar = 1.05457173E-34; //m^2*kg*s^-1
//double h = 0.6711;
double h = 0.7008842;          //relative hubble parameter
double H0;    //s^-1
*/
/*! Partition function for neutrino energy density.
 */
/*double neutrino_partition(double p, void *params)
//double neutrino_partition(double p)
{
  double part;
  double T_neu = *(double *) params;
  double mass_n;
  //extern double T_neu;
  //printf ("T = %f", T_neu);
  mass_n = mass_neu_expan/ktoev;
  part = p*p*sqrt(p*p + mass_n*mass_n)/(exp((p/T_neu) - xi_expan) + 1.0) + p*p*sqrt(p*p + mass_n*mass_n)/(exp((p/T_neu) + xi_expan) + 1.0);
  
  return part;
}
*/
double neutrino_partition(double pot, void *params)
//double neutrino_partition(double p)
{
  double part;
  double T_neu = *(double *) params;
  double mass_n;
  //extern double T_neu;
  //printf ("T = %f", T_neu);
  mass_n = mass_neu_expan/ktoev;
  part = pow(pot, 3)*sqrt(1 + pow(mass_n/(pot*T_neu), 2))/(exp((pot) - xi_expan) + 1.0) + pow(pot, 3)*sqrt(1 + pow(mass_n/(pot*T_neu), 2))/(exp((pot) + xi_expan) + 1.0);
  
  return part;
}


double neutrino_integration(double a)
{
    #define WORKSIZE2 100000
     gsl_function F2;
      gsl_integration_workspace *workspace2;
        
        double unittrans, rrp;
        double inte_result, inte_abserr;
        double roneu;
        double T_neu;

        T_neu = T_neu0/a; 
        H0 = 100. * h * 1000. / ((float)(1E6 * 3.0857E16));
        workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
    unittrans = pow(kb, 4)/((pow(hbar, 3))*(pow(c_light, 3))*2*M_PI*M_PI);
      rrp = 8.0*M_PI*Gr/(cc*3.0*H0*H0);

      F2.function = &neutrino_partition;
      F2.params = &T_neu;
    gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
    //gsl_integration_qag(&F2, 0., 1000000., 0, 1.0e-9, WORKSIZE2, GSL_INTEG_GAUSS61, workspace2, &inte_result, &inte_abserr);
    roneu = inte_result*unittrans*rrp*pow(T_neu, 4);
    gsl_integration_workspace_free(workspace2);
    //printf("h0 = %.20lf, h = %f\n", H0, h);
    return 3*roneu;
}
