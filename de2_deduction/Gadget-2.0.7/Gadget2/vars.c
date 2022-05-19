#include <stdio.h>
#include <gsl/gsl_rng.h>
//#include <fftw3.h>

#include "vars.h"

double frstr_a0;
int Volterra_iteration_num;
int field_reso;

double *a_inte_series;
double *s_inte_series;
double *deltak_inte_series;
double *k1;
double *k2;
double *u;
double a_spacing;
double kcritical;
double neutrino_scheme;
int integration_equation_solver;
int snapno;
int correction_nu_distribution;
int step_index = 1;

double *delta_k_early_real;
double *delta_k_early_im;

