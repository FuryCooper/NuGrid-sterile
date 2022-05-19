#include <stdio.h>
#include <gsl/gsl_rng.h>

#define hbar               1.05457173E-34   //m^2*kg*s^-1
#define kb                 1.3806488E-23    //m^2*kg*s^-2*k^-1
#define Gr                 6.67384E-11      //m^3*kg^-1*s^-2
#define  c                 2.9979e8        //m*s^-1
#define   cc               8.9874e16
#define ktoev               8.617e-5
#define mpc_to_m            3.0857e22
#define ev_to_kg            1.78266191E-36
#define mass_of_sun         1.989E30

//extern double H0;
//extern double rocr;
//extern double h;
//extern double unittrans;
//extern double Omega_m;
//extern double Omega_lambda;
//extern double Omega_nu0;
//extern double rophoton0;


extern double *a_inte_series;
extern double *s_inte_series;
extern double *k1;
extern double *k2;
extern double *u;
extern double a_spacing;
extern int integration_equation_solver;
extern int Volterra_iteration_num;

//extern double kcritical;
extern double neutrino_scheme;

extern int snapno;
extern int step_index;
extern int phi_param;

double neutrino_partition(double pot, void *param);
double neutrino_integration(double a, double m, double xi);
double hubble(double a);
double frstr(double k, double pk0_nu, double pk0_cdm, double pk1_cdm, double *array1, double *array2, double frstr_a0, double frstr_a1, double mass, double xi);
void read_early(void);
double om0(double a0);
double oml(double a0);
double growth_fac_D(double a0, double a1);
double a_to_s(int i, double *array1, double frstr_a0, double frstr_a1);
double perturb_fd_integrand(double pot, void *par);
double perturb_fd_integration(double a, double phi, double m, double xi);
double normal_integrand(double x, void *norm_par);
double cal_xi2(double xi3);
double cal_xi1(double xi2, double xi3);
double numdens_integrand(double x, void *par);
double numdens(double xi);

