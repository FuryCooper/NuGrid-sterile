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

double H0;    //s^-1

struct neu_params {double x; double y; double z;};

double neutrino_partition(double pot, void *par)
{
    struct neu_params * params = (struct neu_params *)par;
    double part;
    double Tneu = (params->x);
    double mass_nu = (params->y);
    double nu_xi = (params->z);
    double mass_n;
    
    mass_n = mass_nu/ktoev;
    part = pow(pot, 3)*sqrt(1 + pow(mass_n/(pot*Tneu), 2))/(exp((pot) - nu_xi) + 1.0) + pow(pot, 3)*sqrt(1 + pow(mass_n/(pot*Tneu), 2))/(exp((pot) + nu_xi) + 1.0);
    //printf("massn %f mass nu %f, ktoev %f part %f Tneu %f\n", mass_n, mass_nu, ktoev, part,Tneu);
    
    return part;
}


double neutrino_integration(double a, double m, double xi)
{
#define WORKSIZE2 100000
    gsl_function F2;
    gsl_integration_workspace *workspace2;
    
    double inte_result, inte_abserr;
    double roneu;
    double Tneu;
    Tneu = Tneu0/a;
    struct neu_params alpha = {Tneu, m, xi};
    
    workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
    F2.function = &neutrino_partition;
    F2.params = &alpha;
    
    gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
    roneu = inte_result * unittrans * pow(Tneu, 4);
    roneu = roneu / rocr;
    gsl_integration_workspace_free(workspace2);

    //return 3*roneu;
    return roneu;
}
