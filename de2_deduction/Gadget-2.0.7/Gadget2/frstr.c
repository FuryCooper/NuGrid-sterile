#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <fftw3.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "vars.h"
#include "allvars.h"

#define relativistic 1

struct my_f_params {double x; double y; double z; double w;};
struct neu_params {double x; double y; double z;};

double phi_integrand(double u, void *par);

double phi(double q, double mass, double xi, double a);
double phi_expansion(double q, double mass, double xi, double a);
double phi_fit(double q, double mass, double a);

double normalize_integrand(double u);

double normalize(void);


double cal_xi2(double xi3){
    double AA, BB, CC, DD, t12, t13, s13, s23, c23, s12, c12, c13, r23;
    double xi2, L3, L2;
    int i;
    
    s12 = sqrt(0.304);
    s23 = sqrt(0.51);
    s13 = sqrt(0.0219);
    c23 = sqrt(1. - 0.51);
    c12 = sqrt(1. - 0.304);
    c13 = sqrt(1. - 0.0219);
    t12 = s12 / c12;
    t13 = s13 / c13;
    
    AA = c23*((1.- t12*t12)*c23 - 2.*s13*s23*t12);
    BB = ((1.- t13*t13)*s23*s23 - t12*t13*t13*c23*(2.*s13*s23+t12*c23));
    CC = s23*((1.-t12*t12)*s23+2.*s13*c23*t12);
    DD = ((1.- t13*t13)*c23*c23 + t12*t13*t13*s23*(2.*s13*c23 - t12*s23));
    
    r23 = (DD - BB) / (AA - CC);
    
    L3 = xi3 * (xi3 * xi3 + M_PI * M_PI);
    L2 = r23 * L3;
    
    xi2 = L2 / (M_PI * M_PI);
    for(i=0;i<10;i++){
        xi2 = L2 / (xi2 * xi2 + M_PI * M_PI);
     //   printf("xi2 %f\n", xi2);
    }
    
    return xi2;
}

double cal_xi1(double xi2, double xi3){
    double s13, s23, c23, s12, c12, c13, r23;
    double xi1, L3, L2, L1;
    int i;
    
    s12 = sqrt(0.304);
    s23 = sqrt(0.51);
    s13 = sqrt(0.0219);
    c23 = sqrt(1. - 0.51);
    c12 = sqrt(1. - 0.304);
    c13 = sqrt(1. - 0.0219);
    
    
    L3 = xi3 * (xi3 * xi3 + M_PI * M_PI);
    L2 = xi2 * (xi2 * xi2 + M_PI * M_PI);
    
    L1 = - (s13*s13*L3 + c13*c13*s12*s12*L2) / (c13*c13*c12*c12);
    xi1 = L1 / (M_PI * M_PI);
    
    for(i=0;i<10;i++){
        xi1 = L1 / (xi1 * xi1 + M_PI * M_PI);
      //  printf("xi1 %f\n", xi1);
    }
    
    return xi1;
}



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
        Tneu = All.Tneu0/a;
        struct neu_params alpha = {Tneu, m, xi};
    
    
        workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
      F2.function = &neutrino_partition;
      F2.params = &alpha;
    
    gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
    roneu = inte_result * All.unittrans * pow(Tneu, 4);
    roneu = roneu / All.rocr;
    
    //printf("mass %f xi %f a %f roneu %f rocr %f unittrans %f inte_result %f\n", m, xi, a, roneu, All.rocr *1e15, All.unittrans*1e15, inte_result);
    
    gsl_integration_workspace_free(workspace2);
    return roneu;
}


double hubble(double a)
{
  double hub, rrp, rom, ror, rolambda, roneu, rototal;

  rom = All.Omega2 * All.rocr / (a * a * a);
  rolambda = All.OmegaLambda * All.rocr;
    if(All.expan_on == 1){
        roneu = neutrino_integration(a, All.mass_1, All.xi_1) + neutrino_integration(a, All.mass_2, All.xi_2) + neutrino_integration(a, All.mass_3, All.xi_3);
    }
    
    if(All.expan_on == 0){
        roneu = neutrino_integration(a, 0., 0.) * 3;
    }
    
  roneu *= All.rocr;
    
  rototal = rom + roneu + rolambda;
  hub = sqrt(8 * M_PI * Gr * rototal / 3) / c;

  return hub;
}


double a_to_s(int i, double *array1, double frstr_a0, double frstr_a1)
{
  int j;
  double temp;
  double a_spacing = (frstr_a1 - frstr_a0) / (double)(All.frstr_interval);
  temp = 0.;
  /*
    if(i == 1){
    temp = 0.5 / pow(array1[i-1], 3) / hubble(array1[i-1]) + 0.5 / pow(array1[i], 3) / hubble(array1[i]);
  }

  else{
    if(i%2 == 0){
      for(j=1;j<i;j++){
            if(j%2 == 0){
                temp = temp + 2. / pow(array1[j], 3) / hubble(array1[j]);
              }
            else{
                temp = temp + 4. / pow(array1[j], 3) / hubble(array1[j]);
            }
          }
      temp = temp + 1. / pow(array1[i], 3) / hubble(array1[i]);
      temp = temp / 3.;
    }

    else{
      for(j=1;j<i-1;j++){
            if(j%2 == 0){
                temp = temp + 2. / pow(array1[j], 3) / hubble(array1[j]);
              }
            else{
                temp = temp + 4. / pow(array1[j], 3) / hubble(array1[j]);
            }
          }
      temp = temp + 1. / pow(array1[i-1], 3) / hubble(array1[i-1]);
      temp = temp / 3.;
      temp = temp + 0.5 / pow(array1[i-1], 3) / hubble(array1[i-1]) + 0.5 / pow(array1[i], 3) / hubble(array1[i]);
    }
  }

   //printf("hubble %f temp %f aspacing %f a %f\n", hubble(array1[i-1])*1e18, temp, a_spacing, array1[i]);
  */
    
    for(j=1;j<=i;j++){
        temp += 0.5 / pow(array1[j-1], 3) / hubble(array1[j-1]) + 0.5 / pow(array1[j], 3) / hubble(array1[j]);
    }
    
   return temp * a_spacing;
}


double frstr(double k, double pk0_nu, double pk0_cdm, double pk1_cdm, double *array1, double *array2, double frstr_a0, double frstr_a1, double mass, double xi)
{
  double *phi2;
  int i, knum, j, m;
  int Volterra_iteration_num = 1;
  double a_spacing = (frstr_a1 - frstr_a0) / (double)(All.frstr_interval);

    double phi0, phi1;
  phi2 = (double*) malloc((All.frstr_interval + 1) * sizeof(double));

  double *k1;
  k1 = (double*) malloc((All.frstr_interval + 1) * sizeof(double));

  double *k2;
  k2 = (double*) malloc((All.frstr_interval + 1) * sizeof(double));
    
    int integration_equation_solver = 1;

  for(i=0;i<All.frstr_interval;i++)
 {
     //printf("arr2 %lf arr1 %lf k %f array 2 %f\n", k*(array2[All.frstr_interval] - array2[i]), array1[i], k, array2[i]);
     if(All.phi_param == 1){
     //phi1[i] = phi(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
         phi2[i] = phi(k*array2[i], mass, xi, array1[i]);
     }
     
     if(All.phi_param == 2){
         //phi1[i] = phi_expansion(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
         phi2[i] = phi_expansion(k*array2[i], mass, xi, array1[i]);
     }
     
     if(All.phi_param == 3){
         //phi1[i] = phi_fit(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
         phi2[i] = phi_fit(k*array2[i], mass, array1[i]);
     }
     
    //printf("phi2 %.10f i %d q %f phi1*Gr*All.rocr%.15f\n", phi2[i], i, array2[i], phi1[i] * Gr * All.rocr);

    k1[i] = 0.;
    k2[i] = 0.;
  }

  //printf("phi loop ended\n");
    if(All.phi_param == 1){
        phi2[All.frstr_interval] = phi(k*array2[All.frstr_interval], mass, xi, array1[All.frstr_interval]);
    }
    
    if(All.phi_param == 2){
        phi2[All.frstr_interval] = phi_expansion(k*array2[All.frstr_interval], mass, xi, array1[All.frstr_interval]);
    }
    
    if(All.phi_param == 3){
        phi2[All.frstr_interval] = phi_fit(k*array2[All.frstr_interval], mass, array1[All.frstr_interval]);
    }

  //calculate 3d delta k
  if(integration_equation_solver == 1){
    //iteration to solve the Volterra equation, for the discrete inte part, trapazoidal rule is temporarily used here
    double f0, f1, m0, m1, f_temp, m_temp, m_temp2, f_temp2;
    
    
    k1[0] = 0.;
    double munit = 4 * M_PI * Gr * All.rocr * a_spacing;
    for(i=1;i<=All.frstr_interval;i++){
        m_temp = 0.;
        for(j=1;j<=i;j++){
        
            double temp_omegam;
        
            temp_omegam = All.Omega0 - All.Omega_nu0_frstr;

            if(All.phi_param == 1){
                phi0 = phi(k*(array2[i] - array2[j-1]), mass, xi, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi(k*(array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(c, 3);
            }
            
            if(All.phi_param == 2){
                phi0 = phi_expansion(k*(array2[i] - array2[j-1]), mass, xi, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi_expansion(k*(array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(c, 3);
                //printf("phi0 %f i %d j %d k %f\n", phi0, i, j, k);
            }
            
            if(All.phi_param == 3){
               phi0 = phi_fit(k*(array2[i] - array2[j-1]), mass, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi_fit(k*(array2[i] - array2[j]), mass, array1[j]) * mpc_to_m / pow(c, 3);
                //printf("phi0 %f i %d j %d ti %f tj %f\n", phi0, i, j, array2[i], array2[j-1]);
            }
            
        m0 = phi0 * array1[j-1] * (array2[i] - array2[j-1]) * temp_omegam * ((pk1_cdm - pk0_cdm) * (j-1)/ All.frstr_interval + pk0_cdm) / (pow(array1[j-1], 3) * hubble(array1[j-1]));
        
        m1 = phi1 * array1[j] * (array2[i] - array2[j]) * temp_omegam * ((pk1_cdm - pk0_cdm) * j / All.frstr_interval + pk0_cdm) / (pow(array1[j], 3) * hubble(array1[j]));
            //printf("phi0 %f phi1 %f\n", phi0, phi1);
         m_temp = m_temp + 0.5 * (m0 + m1);
      }
         m_temp2 = m_temp * munit;
        
      k1[i] = m_temp2 + phi2[i] * pk0_nu;
       // printf("m_temp2 %.10f, phi2[i] %f q %f\n", m_temp2, phi2[i], k*array2[i]);
    }
   //k0 is the higher redshift delta in k space at (kx, ky, kz)

  for(i=0;i<=All.frstr_interval;i++)
  {
    k2[i] = k1[i];
  }

  for(m=0;m<Volterra_iteration_num;m++){
      double onu_temp0, onu_temp1;
    for(i=1;i<=All.frstr_interval;i++){
        f_temp = 0.;
        for(j=1;j<=i;j++){
         
            onu_temp0 = neutrino_integration(array1[j-1], mass, xi);
            onu_temp1 = neutrino_integration(array1[j], mass, xi);
            
            onu_temp0 *= pow(array1[j-1], 3);   //this is because the neutrino_integration func returns Omega in unit of rocr at 0 in terms of percentage
            onu_temp1 *= pow(array1[j], 3);
        
            if(All.phi_param == 1){
                phi0 = phi(k*(array2[i] - array2[j-1]), mass, xi, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi(k*(array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(c, 3);
            }
            
            if(All.phi_param == 2){
                phi0 = phi_expansion(k*(array2[i] - array2[j-1]), mass, xi, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi_expansion(k*(array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(c, 3);
            }
            
            if(All.phi_param == 3){
                phi0 = phi_fit(k*(array2[i] - array2[j-1]), mass, array1[j-1]) * mpc_to_m / pow(c, 3);
                phi1 = phi_fit(k*(array2[i] - array2[j]), mass, array1[j]) * mpc_to_m / pow(c, 3);
            }
            
            
            f0 = phi0 * array1[j-1] * (array2[i] - array2[j-1]) * onu_temp0 * k2[j-1] / (pow(array1[j-1], 3) * hubble(array1[j-1]));
            f1 = phi1 * array1[j] * (array2[i] - array2[j]) * onu_temp1 * k2[j] / (pow(array1[j], 3) * hubble(array1[j]));
      
            f_temp = f_temp + 0.5 * (f0 + f1);
        }
        f_temp2 = f_temp * munit;

      k2[i] = f_temp2 + k1[i];
      //printf("ftemp %lf k2%lf k1%lf f0%lf\n", f_temp2, k2[i], k1[i], f0);
     }
    //printf("Volterra eqn no %d, delta k nv %.8f k1 %.8f knv0 %lf k1part %.8f\n", m, k2[All.frstr_interval], k1[All.frstr_interval], pk0_nu, phi2[All.frstr_interval] * pk0_nu);
    }
    
  }
  return k2[All.frstr_interval];
}


double phi_integrand(double u, void *par)
{
    double up, down1, down2, T;
    
    int unitless_inte = 1;
    struct my_f_params * params = (struct my_f_params *)par;
    
    double a = (params->x);
    double q = (params->y);
    double mass = (params->z);
    double xi = (params->w);
    
    T = All.Tneu0/a;
    double con_ul = T * ktoev / mass;
    
    if(unitless_inte == 1){
        
        if(relativistic == 0){
            up = u * (sin(u * q * con_ul) + 1.);
            down1 = exp(u - xi) + 1.;
            down2 = exp(u + xi) + 1.;
            
        }
        if(relativistic == 1){
            up = u / sqrt(1. - u*u*con_ul*con_ul) * (sin(u / sqrt(1. - u*u*con_ul*con_ul) * q * con_ul) + 1.);
            down1 = exp(u / sqrt(1. - u*u*con_ul*con_ul) - xi) + 1.;
            down2 = exp(u / sqrt(1. - u*u*con_ul*con_ul) + xi) + 1.;
        }
        //up = u * sin(u * q * con_ul);
        //expand in taylor
        
    }
    
    return up / down1 + up / down2;
    
}


double phi_to_deduct(double u, void *par)
{
    double up, down1, down2, T;
    
    struct my_f_params * params = (struct my_f_params *)par;
    
    double a = (params->x);
    double q = (params->y);
    double mass = (params->z);
    double xi = (params->w);
    
    T = All.Tneu0/a;
    double con_ul = T * ktoev / mass;
    //up = u * sin(u * q * con_ul);
    if(relativistic == 0){
        up = u;
        down1 = exp(u - xi) + 1.;
        down2 = exp(u + xi) + 1.;
    }
    if(relativistic == 1){
        up = u / sqrt(1. -  u*u*con_ul*con_ul);
        down1 = exp(u / sqrt(1. - u*u*con_ul*con_ul) - xi) + 1.;
        down2 = exp(u / sqrt(1. - u*u*con_ul*con_ul) + xi) + 1.;
    }
    
    return up / down1 + up / down2;
}


/*double normal_integrand(double x, void *norm_par){
    
    double part1, part2;
    double Tneu = *(double *) norm_par;
    double con_ul = Tneu * ktoev / All.mass_nu_frstr;
    
    if(relativistic == 0){
        part1 = x * x / (exp(x - xi) + 1.);
        part2 = x * x / (exp(x + xi) + 1.);
    }
    if(relativistic == 1){
        x = x / sqrt(1. - x*x*con_ul*con_ul);
        part1 = x * x / (exp(x - xi) + 1.);
        part2 = x * x / (exp(x + xi) + 1.);
    }
    return (part1 + part2);
    
}*/

double phi(double q, double mass, double xi, double a)
{
  double inte_result, inte_abserr;
  double inte_result2, inte_abserr2;
    double normal_factor;
    double TT = All.Tneu0 / a;
#define WORKSIZE2 100000
     gsl_function F2;
     gsl_function F3;
    gsl_function F4;
      gsl_integration_workspace *workspace2;
      gsl_integration_workspace *workspace3;
    gsl_integration_workspace *workspace4;
       
       struct my_f_params alpha = {a, q, mass, xi};

  workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
  workspace3 = gsl_integration_workspace_alloc(WORKSIZE2);
workspace4 = gsl_integration_workspace_alloc(WORKSIZE2);
  F2.function = &phi_integrand;
  F2.params = &alpha;

  F3.function = &phi_to_deduct;
  F3.params = &alpha;
    
    //F4.function = &normal_integrand;
    //F4.params = &TT;
    
  double con_ul = (All.Tneu0 / a) * ktoev / mass;
  //gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-2, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
  gsl_integration_qag(&F2, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace2, &inte_result, &inte_abserr);
  gsl_integration_qag(&F3, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace3, &inte_result2, &inte_abserr2);
    //gsl_integration_qag(&F4, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace4, &normal_factor, &inte_abserr2);
  //gsl_integration_qagiu(&F3, 0.0, 0, 1.0e-2, WORKSIZE2, workspace3, &inte_result2, &inte_abserr2);
  gsl_integration_workspace_free(workspace2);
  gsl_integration_workspace_free(workspace3);
    gsl_integration_workspace_free(workspace4);
  //printf("phi check point 2 %.20f %lf\n", inte_result, q);
    
    int unitless_inte = 1;
  inte_result = inte_result - inte_result2;
  double zeta3 = 1.202056;

  if(unitless_inte == 1){
    //inte_result = inte_result * 4. / 3. / q / con_ul / zeta3;
      inte_result = inte_result / q / con_ul / normal_factor;
    }
    //question why integrate to 2.?
  return inte_result;

}


double phi_expansion(double q, double mass, double xi, double a){
    double Tneu = All.Tneu0 / a;
    double zeta3 = 1.202056;
    int i, n;
    double sum, temp, x, sum2;
    double con_ul = Tneu * ktoev / mass;
    
    sum = 0.;
    sum2 = 0.;
    x = q * con_ul;
    n = 500;
    
    if(q == 0.){
        sum2 == 0.;
    }
    
    else{
        
    if(xi == 0.){
        for(i=1;i<=n;i++){
            temp = pow(-1., i+1) * i / pow((i*i + x*x), 2);
            sum += temp;
        }
        sum2 = sum;
        sum2 *= 4. / (3. * zeta3);
    }

    if(xi != 0.){
        
        double temp_dens;
        
        if(xi == All.xi_1){
            temp_dens = All.numdens1;
        }
        if(xi == All.xi_2){
            temp_dens = All.numdens2;
        }
        if(xi == All.xi_3){
            temp_dens = All.numdens3;
        }
        
        xi = abs(xi);
        double c1, a1, a2, a3, a4, sum, a5, temp1, temp2, err;
        a1 = 0.;
        a2 = 0.;
        a3 = 0.;
        a4 = 0.;
        a5 = 0.;
        
        for(i=1;i<=n;i++){
            a1 += pow(-1., i+1) * i / (i*i + x*x);
            a2 += pow(-1., i+1) * x / (i*i + x*x);
            a3 += pow(-1., i+1) * (i*i - x*x) / pow((i*i + x*x), 2);
            a4 += pow(-1., i+1) * (2.*x*i) / pow((i*i + x*x), 2);
            a5 += pow(-1., i+1) * 2. * i * x * exp(-i*xi) / pow((i*i + x*x), 2);
        }
        
        struct my_f_params alpha = {a, q, mass, xi};
        #define WORKSIZE2 100000
        gsl_function F2;
        gsl_function F3;
        
        gsl_integration_workspace *workspace2;
        gsl_integration_workspace *workspace3;
        workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
        workspace3 = gsl_integration_workspace_alloc(WORKSIZE2);
        
        F2.function = &phi_integrand;
        F2.params = &alpha;
        
        F3.function = &phi_to_deduct;
        F3.params = &alpha;

        
        gsl_integration_qag(&F2, 0.0, xi, 0, 1.0e-6, WORKSIZE2, GSL_INTEG_GAUSS61, workspace2, &temp1, &err);
        gsl_integration_qag(&F3, 0.0, xi, 0, 1.0e-6, WORKSIZE2, GSL_INTEG_GAUSS61, workspace3, &temp2, &err);
        gsl_integration_workspace_free(workspace2);
        gsl_integration_workspace_free(workspace3);
        
        c1 = temp1 - temp2;
        
        /*int num_inte_step = 30;
        c1 = 0.;
        //double *num_inte_temp;
        double temp_x0, temp_x1, temp_func;
        //num_inte_temp = (double*) malloc((num_inte_step + 1) * sizeof(double));
        
        for(i=0;i<num_inte_step;i++){
            temp_x0 = xi * i / num_inte_step;
            temp_x1 = xi * (i+1) / num_inte_step;
            temp_func = 0.25 * (temp_x0 * sin(x * temp_x0) / (exp(temp_x0 - xi) + 1) + temp_x1 * sin(x * temp_x1) / (exp(temp_x1 - xi) + 1));
            temp_func += 0.25 * (temp_x0 * sin(x * temp_x0) / (exp(temp_x0 + xi) + 1) + temp_x1 * sin(x * temp_x1) / (exp(temp_x1 + xi) + 1));
            c1 += temp_func * xi / num_inte_step;
        }*/
        
        sum = c1 + xi * cos(x*xi) * a2 + xi * sin(x*xi) * a1 + cos(x*xi) * a4 + sin(x*xi) *a3 + a5;
        //printf("c1 %f a2 %f a1 %f a3 %f a4 %f a5 %f\n", c1 * 1e10, xi * cos(x*xi)*a2*1e10, xi * sin(x*xi)*a1*1e10, sin(x*xi)*a3*1e10, cos(x*xi)*a4*1e10, a5*1e10);
        sum2 = sum / 2. / x / 2.; //the second 2 comes from the nu-antinu pair, that in the normalization part cancelled
        //printf("sum %f x %f\n", sum, x);
         sum2 *= 4. / temp_dens;
        
    }
       
    }
    

    //sum2 *= 4. / (3. * zeta3) * All.numdens0 / All.numdens1;
    
    
    if(sum2 > 1.){
        sum2 = 0.99999;
    }

    return sum2;
    
}

double phi_fit(double q, double mass, double a){
    double Tneu = All.Tneu0 / a;
    double zeta3 = 1.202056;
    double con_ul = Tneu * ktoev / mass;
    double fit, x;
    
    x = q * con_ul;

    fit = (1. + 0.0168*x*x + 0.0407*pow(x, 4)) / (1. + 2.1734*x*x + 1.6787*pow(x, 4.1811) + 0.1467*pow(x, 8));
    
    return fit;
    
}

double numdens_integrand(double x, void *par){
    
    double xi = *(double *) par;
    return x * x /(exp(x - xi) + 1.) + x * x /(exp(x + xi) + 1.);
    
}


double numdens(double xi){
#define WORKSIZE2 100000
    gsl_function F;
    gsl_integration_workspace *workspace2;
    
    double inte_result, inte_abserr;
    
    workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
    F.function = &numdens_integrand;
    F.params = &xi;
    
    gsl_integration_qagiu(&F, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
    gsl_integration_workspace_free(workspace2);
    
    return inte_result;
    
}

/*
double oml(double a0)
{
    double om0a0;
    
    om0a0 = All.Omega0 * pow(a0, -3.) / (All.Omega0 * pow(a0, -3.) + All.OmegaLambda + neutrino_integration(a0, All.mass_nu_expan, All.xi_expan));
    //printf("a0 %f, oml %f, Omegam %f, roneu %f\n", a0, om0a0, All.Omega0, neutrino_integration(a0));
    return om0a0;
}



double om0(double a0)
{
    double om0l0;
    
    om0l0 = All.OmegaLambda / (All.Omega0 * pow(a0, -3.) + All.OmegaLambda + neutrino_integration(a0, All.mass_nu_expan, All.xi_expan));
    
    return om0l0;
}


double growth_fac_D(double a0, double a1)
{
    double D0, D1;
    
    D0 = a0 * 2.5 * om0(a0) / ( pow(om0(a0), 4./7.) - oml(a0) + (1 + om0(a0) / 2.) * (1 + oml(a0) / 70.) );
    D1 = a1 * 2.5 * om0(a1) / ( pow(om0(a1), 4./7.) - oml(a1) + (1 + om0(a1) / 2.) * (1 + oml(a1) / 70.) );
    //printf("om0 %f oml %f D0 %f a0 %f\n", om0(a0), oml(a0), D0, a0);
    return D1 / D0;
}
 */
