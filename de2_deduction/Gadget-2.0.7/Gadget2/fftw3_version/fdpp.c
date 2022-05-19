#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <fftw3.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
//#include        <rfftw_mpi.h>
#include     <drfftw_mpi.h>	/* double precision FFTW */
#include     <srfftw_mpi.h>

#include "vars.h"
#include "allvars.h"


//#define field_resolution 256
//#define All.BoxSizeize 200000
#define M_PI 3.14159265358979323846

struct fdpp_params {double x; double y;}; 

double perturb_fd_integrand(double pot, void *par){
   

    double part;
    double mass_n;

    struct fdpp_params * params = (struct fdpp_params *)par;
    double Tneu = (params->x);
    double phi = (params->y);

    mass_n = mass_neu/ktoev;
    part = pow(pot, 2) / (exp((pot) - nu_xi + phi * mass_n / Tneu / cc) + 1.0) + pow(pot, 2) / (exp((pot) + nu_xi + phi * mass_n / Tneu / cc) + 1.0);
    //part = pow(pot, 3) * sqrt(1 + pow(mass_n/(pot*Tneu), 2)) / (exp((pot) - nu_xi + phi * mass_n / Tneu / cc) + 1.0) + pow(pot, 3) * sqrt(1 + pow(mass_n/(pot*Tneu), 2)) / (exp((pot) + nu_xi + phi * mass_n / Tneu / cc) + 1.0);
    
    return part;

}



double perturb_fd_integration(double a, double phi)
{
    double Tneu = Tneu0/a;
    struct fdpp_params alpha = {Tneu, phi};

    #define WORKSIZE2 100000
     gsl_function F2;
      gsl_integration_workspace *workspace2;
        
        double inte_result, inte_abserr;
        double roneu;

        workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
      
      F2.function = &perturb_fd_integrand;
      F2.params = &alpha;
    gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
    roneu = inte_result*unittrans*pow(Tneu, 3) / kb;
    //roneu = inte_result*unittrans*pow(Tneu, 4);
    gsl_integration_workspace_free(workspace2);

    return 3*roneu;
}



void phi_perturb_correction(fftw_real *workspace){
    
    int xi, yi, zi, xtemp, ytemp, ztemp;
    double k, kx, ky, kz;
    double average;
    double Omega_nu0 = 0.007396;
    average = (double)(NumPart) / (double)(pow(PMGRID+1, 3));
    
    int nmid = PMGRID / 2;
    
    fftw_complex *fluc_density_field;
    fluc_density_field = (fftw_complex*) fftw_malloc((PMGRID+1)*(PMGRID+1)*(PMGRID+1) * sizeof(fftw_complex));
    fftw_complex *out;
    out = (fftw_complex*) fftw_malloc((PMGRID+1)*(PMGRID+1)*(PMGRID+1) * sizeof(fftw_complex));
    
    fftw_plan fluc_density_field_plan;
    fluc_density_field_plan = fftw_plan_dft_3d((PMGRID+1), (PMGRID+1), (PMGRID+1), fluc_density_field, out, 1, FFTW_MEASURE);

    fftw_complex *delta_ro;
    delta_ro = (fftw_complex*) fftw_malloc((PMGRID+1)*(PMGRID+1)*(PMGRID+1) * sizeof(fftw_complex));


        for(xi=0;xi<=PMGRID;xi++){
            for(yi=0;yi<=PMGRID;yi++){
                for(zi=0;zi<=PMGRID;zi++){
                    //fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].re = (workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)] - average)/average;
                    fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].re = (workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)])/average;
                    fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].im = 0.;
                    /*if(fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].re < 0.){
                        
                        printf("dens = %f fluc = %f ave = %f xyz = %d %d %d\n", workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)], fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)], average, xi, yi, zi);
                    }*/
                }
            }
        }
        
        fftw_execute(fluc_density_field_plan);
        
        fftw_complex *phik;
        fftw_complex *phi_real;
        
        phik = (fftw_complex*) malloc((PMGRID + 1)*(PMGRID + 1)*(PMGRID + 1) * sizeof(fftw_complex));
        phi_real = (fftw_complex*) malloc((PMGRID + 1)*(PMGRID + 1)*(PMGRID + 1) * sizeof(fftw_complex));
        
        fftw_plan plan_of_phi;
        plan_of_phi = fftw_plan_dft_3d((PMGRID+1), (PMGRID+1), (PMGRID+1), phik, phi_real, -1, FFTW_MEASURE);
        
        for(xi=0;xi<=PMGRID;xi++){
            for(yi=0;yi<=PMGRID;yi++){
                for(zi=0;zi<=PMGRID;zi++){
                    
                    phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = 0.;
                    phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = 0.;
                    phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = 0.;
                    phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = 0.;
                    
                }
            }
        }
        
        
        
        for(xi=0;xi<=PMGRID;xi++){
            for(yi=0;yi<=PMGRID;yi++){
                for(zi=0;zi<=PMGRID;zi++){
                    
                    if(xi>nmid){
                        xtemp = PMGRID + 1 - xi;
                    }
                    else{
                        xtemp =  - xi;
                    }
                    
                    if(yi>nmid){
                        ytemp = PMGRID + 1 - yi;
                    }
                    else{
                        ytemp =  - yi;
                    }
                    
                    if(zi>nmid){
                        ztemp = PMGRID + 1 - zi;
                    }
                    else{
                        ztemp =  - zi;
                    }
                    
                    kx = (2*M_PI*1000/All.BoxSize) * xtemp;  //MPc ^ -1
                    ky = (2*M_PI*1000/All.BoxSize) * ytemp;
                    kz = (2*M_PI*1000/All.BoxSize) * ztemp;
                    
                    k = sqrt(kx*kx + ky*ky + kz*kz);
                    
                    if(k == 0){
                        phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = 0.;
                        phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = 0.;
                    }
                    
                    else{
                        phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = (0. - out[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re) * rocr / cc * All.Omega0 * 4. * M_PI * All.G / ((All.Time * All.Time * All.Time) * (k * k / mpc_to_m / mpc_to_m) * pow((PMGRID+1), 3));
                        phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = (0. - out[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im) * rocr / cc * All.Omega0 * 4. * M_PI * All.G / ((All.Time * All.Time * All.Time) * (k * k / mpc_to_m / mpc_to_m) * pow((PMGRID+1), 3));
                        //printf("phi = %lf, k = %lf out = %lf norm %f\n", phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re, mpc_to_m * mpc_to_m / k / k / 1e40, out[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re * rocr * All.Omega0 * 4. * M_PI * All.G * 1e15, pow((PMGRID+1), 3));
                    }
                }
            }
        }
        
        
        fftw_execute(plan_of_phi);
        
        /*
        for(xi=0;xi<=PMGRID;xi++){
            for(yi=0;yi<=PMGRID;yi++){
                for(zi=0;zi<=PMGRID;zi++){
                    
                    //phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re;
                    //phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im;
                    //phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im = phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im / pow((PMGRID+1), 3);
                    //printf("phi_real imaginary %f real %f phik %lf\n", phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].im, phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re, phik[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re);
                    //printf("phi/T %f phi %f\n", phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re * 0.1 / ktoev / cc / 1.7, phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re);
                }
            }
        }*/
        
        double ev_to_kg = 1.78266191E-36;
        double mass_of_sun = 1.989E30;
        double average2;
        double rosum = 0.;
        average2 = rocr * All.Omega0 * pow((200 * mpc_to_m / h), 3) / cc / mass_of_sun / 1e10 * h / 33.5;
        
        for(xi=0;xi<=PMGRID;xi++){
            for(yi=0;yi<=PMGRID;yi++){
                for(zi=0;zi<=PMGRID;zi++){
                    
                    delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = mass_neu * ev_to_kg * perturb_fd_integration(All.Time, phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re) * pow((mpc_to_m * All.BoxSize / 1000. / (PMGRID+1)), 3) / h / h / mass_of_sun / 1E10 / P[1].Mass;
                    
                    rosum = rosum + delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re / average;
                    //delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = mass_neu * ev_to_kg * perturb_fd_integration(All.Time, 0.) * pow((mpc_to_m * All.BoxSize / 1000. / PMGRID), 3) * h * h / mass_of_sun / 1E10 / P[1].Mass;
                    //delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re = perturb_fd_integration(All.Time, 0.) * pow((mpc_to_m * All.BoxSize / 1000. / PMGRID), 3) * h * h / mass_of_sun / 1E10 / P[1].Mass / cc;
                    //double cdm = workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)];
                    //double neuneu = mass_neu * ev_to_kg * perturb_fd_integration(All.Time, phi_real[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re) * cc;
                    
                    /*if(cdm > 1){
                     printf("delta_ro = %lf, cdm mass density %lf mass %lf average %lf ro %lf omega ratio %lf\n", delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re / average, cdm, P[1].Mass, average, Omega_nu0);
                     }*/
                    
                }
            }
        }
        rosum = rosum / pow(PMGRID+1, 3);
        printf("rosum/average %lf Omega_nu0 / All.Omega0 %lf\n", rosum, neutrino_integration(1.) / rocr / All.Omega0);
    
    int correction_nu_distribution = 2;
    for(xi=0;xi<=PMGRID;xi++){
        for(yi=0;yi<=PMGRID;yi++){
            for(zi=0;zi<=PMGRID;zi++){
                double roro = delta_ro[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)].re;
                if(correction_nu_distribution == 1){
                    //fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].re = (density_field[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)] * All.Omega0 / (All.Omega0 + Omega_nu0) + (roro + 1.) * average * Omega_nu0 / (All.Omega0 + Omega_nu0) - average)/average;
                    workspace[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)] = workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)] + (roro + 1.) * average * Omega_nu0 / All.Omega0;
                }
                
                if(correction_nu_distribution == 2){
                    workspace[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)] = workspace[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)] * All.Omega0 / (All.Omega0 + Omega_nu0) + roro;
                }
                //if(fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)] != 0.){
                
                //printf("dens = %f fluc = %f ro %lf ave = %f xyz = %d %d %d\n", density_field[zi+(PMGRID+1)*(yi + (PMGRID+1)*xi)], fluc_density_field[zi + (PMGRID+1)*(yi + (PMGRID+1)*xi)].re, roro, average, xi, yi, zi);
                //}
            }

        }
}
}
