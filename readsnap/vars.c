#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>

#include "vars.h"

int ThisTask;
int snapno;
char output[200];
char output1[200];
char output2[200];
char output3[200];
char input_path[200];
char nu_txt[200];
double fnu;
