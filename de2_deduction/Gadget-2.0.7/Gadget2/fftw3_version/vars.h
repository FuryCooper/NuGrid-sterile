#include <stdio.h>
#include <gsl/gsl_rng.h>
//#include <fftw3.h>

extern double nu_xi; 						//neutrino chemical potential over temperature
extern double Tneu0;					//neutrino tempearture at z=0
extern double mass_neu;				//neutrino mass
extern double T_cmb;	

extern double hbar;
extern double kb;
extern double Gr;
extern double c;
extern double cc;
extern double ktoev;

extern double H0;
extern double rocr;
extern double h;
extern double roneu0;
extern double mnu_frstr;
extern double rrp;
extern double unittrans;
extern double Omega_m;
extern double Omega_lambda;
extern double Omega_nu0;
extern double rophoton0;
extern double mpc_to_m;
extern double normal_factor;

extern double frstr_a0;
extern int frstr_interval;

extern double *a_inte_series;
extern double *s_inte_series;
extern double *deltak_inte_series;
extern double *k1;
extern double *k2;
extern double *u;
extern double *initial_pk;
extern double *inte_matrix;
extern double a_spacing;
extern int integration_equation_solver;
extern int Volterra_iteration_num;

extern int Nbins;
extern double *delta_k_early_real;
extern double *delta_k_early_im;
extern double kcritical;

extern int field_reso;
extern int ThisTask;
extern int snapno;
extern int correction_nu_distribution;

extern char addfload[100];
extern char txtname[100];
extern char f_load_early[200];
extern char f_load_early2[200];

double neutrino_partition(double pot, void *param);
double neutrino_integration(double a);
double hubble(double a);
double frstr(double k, double k0, double *array1, double *array2);
void read_early(void);
double om0(double a);
double oml(double a);
double growth_fac_D(double a);
double a_to_s(int i, double *array1);
double perturb_fd_integrand(double pot, void *par);
double perturb_fd_integration(double a, double phi);
