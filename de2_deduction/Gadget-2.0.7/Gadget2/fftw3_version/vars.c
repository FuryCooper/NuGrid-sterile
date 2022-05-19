#include <stdio.h>
#include <gsl/gsl_rng.h>
//#include <fftw3.h>

#include "vars.h"

double nu_xi; 						//neutrino chemical potential over temperature
double Tneu0;					//neutrino tempearture at z=0
double mass_neu;				//neutrino mass
double T_cmb;	
double mnu_frstr;
double Omega_m;
double Omega_lambda;
double Omega_nu0;

double frstr_a0;
int frstr_interval;
int Volterra_iteration_num;
int field_reso;

double hbar;
double kb;
double Gr;
double c;
double cc;
double ktoev;
double H0;
double rocr;
double h;
double roneu0;
double rrp;
double unittrans;
double rophoton0;
double mpc_to_m;
double normal_factor;

double *a_inte_series;
double *s_inte_series;
double *deltak_inte_series;
double *k1;
double *k2;
double *u;
double *initial_pk;
double *inte_matrix;
double a_spacing;
double kcritical;
int integration_equation_solver;
int snapno;
int correction_nu_distribution;

char addfload[100];
char txtname[100];

char f_load_early[200];
char f_load_early2[200];

int Nbins;
double *delta_k_early_real;
double *delta_k_early_im;

int field_reso;
int ThisTask;
