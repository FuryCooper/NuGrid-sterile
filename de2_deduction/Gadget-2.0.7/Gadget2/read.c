#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "vars.h"

int particle_number, one_d_num, one_d_box;
int box_s;
#define M_PI 3.14159265358979323846
#define Nbins 85
#define field_resolution 256
//#define box_size 40000
#define box_size 200000
#define kbininterval 1.05
#define start_k 0.05
#define nnn 257*257*257
//#define kcritical 0.8
//#define snapno 6
#define sreso 100
//#define correction_nu_distribution 2
//char addfload[] = "/home/zzc/Gadget_output/jl13/standard/123460";
//char txtname[] = "rela_z16";
//char addfload = "/home/zzc/Gadget_output/jl13/standard/123460/z25"

int fractal_density = 1;
int reading_today = 1;

double pk_nor;      //normalization constant for power spectrum
double kn = M_PI / (box_size / field_resolution);

struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header1;



int NumPart, Ngas;

struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
} *P;

int *Id;

double Time, Redshift;



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  printf("checkpoint0\n");

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  //MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    if(argc < 2)
    {
      if(ThisTask == 0)
 {
    fprintf(stdout, "\nParameters are missing.\n");
    fprintf(stdout, "Call with <ParameterFile>\n\n");
  }
      MPI_Finalize();
      exit(0);
    }

  read_parameterfile(argv[1]);
  
  unittrans = pow(kb, 4)/((pow(hbar, 3))*(pow(c, 3))*2*M_PI*M_PI);

  printf("checkpoint0.5\n");
  H0 = 100. * h * 1000. / ((float)(1E6 * 3.0857E16));
  rrp = 8.0*M_PI*Gr/(cc*3.0*H0*H0);
  rocr = 1. / rrp ;
  rophoton0 = M_PI * M_PI * pow(T_cmb, 4) * pow(kb, 4) / (15. * pow(c, 3) * pow(hbar, 3));
  //printf("checkpoint0.6 %f %f %f %f sdsdsad%f\n", rocr*1e18, unittrans, M_PI, Gr, H0);
  printf("rocr %f rrp %f\n", rocr*1e8, rrp);
  printf("nu_xi %f\n", nu_xi);

  //roneu0 = neutrino_integration(1.);
  //printf("checkpoint0.7 roneu0 %lf\n", roneu0/rocr);
  roneu0 = rocr * 0.007396 * mnu_frstr / 0.1;
  printf("checkpoint0.7 roneu0 %lf\n", roneu0/rocr);
  Omega_nu0 = roneu0 / rocr;
  Omega_lambda = Omega_lambda - Omega_nu0;

  char path[200], input_fname[200], basename[200];
  int type, snapshot_number, files, kk;
  int start_time, finish_time;
    box_s = box_size;
    snapshot_number = snapno;		/* number of snapshot */
    files = 1;			/* number of files per snapshot */
    printf("checkpoint1\n");
   
    if (reading_today == 1){
        sprintf(path, addfload);
         sprintf(basename, "snapshot");
        sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
    }
    else {sprintf(path, "/home/zzc/ICs/mar1");
        sprintf(basename, "123457-random-ics");
        sprintf(input_fname, "%s/%s", path, basename);
    }
     printf("checkpoint2\n");
  
  start_time = clock();
  load_snapshot(input_fname, files);

  printf("particle number is %d mass = %f\n", particle_number, P[1].Mass);
  pk_nor = pow(box_size/1000, 3)/(pow(field_resolution, 6));

  one_d_num = (int)(pow(particle_number, 1./3.)+0.5);
  one_d_box = field_resolution;
  printf("one_d_box = %d\n", one_d_box);
    

  //set the frstr correction matrix part
  double *a_inte_series;
  a_inte_series = (double*) malloc((frstr_interval + 1) * sizeof(double));

  double *s_inte_series;
  s_inte_series = (double*) malloc((frstr_interval + 1) * sizeof(double));

  double *deltak_inte_series;
  deltak_inte_series = (double*) malloc((frstr_interval + 1) * sizeof(double));

  double *u;
  u = (double*) malloc((frstr_interval + 1) * sizeof(double));

  printf("checkpoint3\n");

  int i;
  double a_spacing = (1. - frstr_a0) / frstr_interval;

  for(i=0;i<=frstr_interval;i++){
    a_inte_series[i] = frstr_a0 + i * a_spacing;
  }

  printf("checkpoint3.3 %lf\n", a_spacing);

    s_inte_series[0] = 0.;
    printf("checkpoint3.3.1 %d\n", frstr_interval);

    //normal_factor = normalize();

    //printf("%f %f\n", hubble(0.001)*1e15, 0.001);
/*
    for(i=0;i<=frstr_interval;i++){
    printf("%f, %f\n", a_inte_series[i], hubble(a_inte_series[i])*1e10);
    
  }*/

    for(i=1;i<=frstr_interval;i++){
    s_inte_series[i] = a_to_s(i, a_inte_series) * c / mpc_to_m;
    //printf("a %lf s %lf\n", a_inte_series[i], s_inte_series[i]);
  }

  s_inte_series[0] = 0. ;//2 * s_inte_series[1] - s_inte_series[2];
  //printf("s0 %f %f %f\n", s_inte_series[0], c, mpc_to_m);

  printf("checkpoint3.5\n");

  reordering();			/* call this routine only if your ID's are set properly */
  unit_conversion();		/* optional stuff */
  do_what_you_want(a_inte_series, s_inte_series);
  finish_time = clock();
  printf("CLOCKS_PER_SEC = %d\n", CLOCKS_PER_SEC);
  printf("duration = %f\n", (double)((finish_time - start_time)/CLOCKS_PER_SEC));

}


double maxoftwo(double a, double b){
    double max;
    if(a<b){
        max = b;
    }
    else{
        max = a;
    }
    return max;
}

double w(double k1, double k2, double k3){
  double up, down;
  up = sin(M_PI * k1 / (2. * kn)) * sin(M_PI * k2 / (2. * kn)) * sin(M_PI * k3 / (2. * kn));
  down = (M_PI * k1 / (2. * kn)) * (M_PI * k2 / (2. * kn)) * (M_PI * k3 / (2. * kn));
  //printf("up %f down %f k1 %f kn %f k1/kn %f \n", sin(M_PI * k1 / (2. * kn)), (M_PI * k2 / (2. * kn)), k1, kn, k1/kn);
  return (up/down) * (up/down);
}



/* here the particle data is at your disposal 
 */

int do_what_you_want(double *array_a, double *array_s)
{
    int i, xi, yi, zi, b, j, iteration_count;
    int nx = one_d_box + 1;
    int ny = one_d_box + 1;
    int nz = one_d_box + 1;
	double l01, l02, l11, l12, l21, l22, unit_weight, max_weight, box_unit_length, average, k, kx, ky, kz;
    

    int xtemp, ytemp, ztemp;
    int nmid = one_d_box/2;

    fftw_complex *fluc_density_field;
    double *density_field;
    double a_weight[8];
    fftw_complex *out;
    density_field = (double*) malloc((one_d_box+1)*(one_d_box+1)*(one_d_box+1) * sizeof(double));

    double *kbins;
    kbins = (double*) malloc(Nbins * sizeof(double));

    double kbins0[Nbins];
    double kbins1[Nbins];
    double kbinsgenerate[Nbins+1];
    int count[Nbins];
    int count1[Nbins];
    double kerr[Nbins];
    double pk[Nbins], pk_corrected[Nbins], kerr_corrected[Nbins];

    double *load_early_re;
    load_early_re = (double*) malloc(nnn * sizeof(double));

    double *load_early_im;
    load_early_im = (double*) malloc(nnn * sizeof(double)); 

    fftw_complex *load_early;
    load_early = (fftw_complex*) fftw_malloc(nnn * sizeof(fftw_complex));
    
    fluc_density_field = (fftw_complex*) fftw_malloc((one_d_box+1)*(one_d_box+1)*(one_d_box+1) * sizeof(fftw_complex));
    out = (fftw_complex*) fftw_malloc(nx*ny*nz * sizeof(fftw_complex));

    fftw_complex *delta_ka1;
    delta_ka1 = (fftw_complex*) fftw_malloc((one_d_box+1)*(one_d_box+1)*(one_d_box+1) * sizeof(fftw_complex));

    fftw_complex *delta_ro;
    delta_ro = (fftw_complex*) fftw_malloc((one_d_box+1)*(one_d_box+1)*(one_d_box+1) * sizeof(fftw_complex));

    fftw_plan plan_of_nu;
    plan_of_nu = fftw_plan_dft_3d((one_d_box+1), (one_d_box+1), (one_d_box+1), delta_ka1, delta_ro, -1, FFTW_MEASURE);

    fftw_plan fluc_density_field_plan;
    fluc_density_field_plan = fftw_plan_dft_3d((one_d_box+1), (one_d_box+1), (one_d_box+1), fluc_density_field, out, 1, FFTW_MEASURE);

    sprintf(f_load_early, "%s/k_space_3d_re.txt", addfload);
    sprintf(f_load_early2, "%s/k_space_3d_im.txt", addfload);

    printf("in what you want check 1\n");
    kbinsgenerate[0] = start_k;
    
    for(i=1;i<Nbins+1;i++){
        //kbinsgenerate[i] = kbinsgenerate[0]+ kbininterval*i;
        kbinsgenerate[i] = kbinsgenerate[i-1] * kbininterval;
    }
    /*initialize bins arrays*/
    for(i=0;i<Nbins;i++){
        kbins0[i] = kbinsgenerate[i];
        kbins1[i] = kbinsgenerate[i+1];
        count[i] = 0;
        count1[i] = 0;
        kerr[i] = 0.;
        pk[i] = 0.;
        
    }
    
    for(i=0;i<Nbins;i++){
        //kbins[i] = (kbins0[i]+kbins1[i])/2.;
        kbins[i] = sqrt(kbins0[i]*kbins1[i]);
    }
    
    /*test the size of snapshot*/
    double temp_posi;
    double temp2;
    for(i=1;i<=particle_number;i++){
        
        temp2 = maxoftwo(temp_posi, P[i].Pos[0]);
        temp_posi = temp2;
    }

    printf("kbinsmin = %f kbinsmax = %f\n", kbins[0], kbins[Nbins-1]);

    box_unit_length = (double)(box_s)/(double)(one_d_box);
    printf("boxunitlength = %f\n", box_unit_length);
    /*initialize density field array*/

    for(i=0;i<=7;i++){
        a_weight[i] = 0.;
    }

    for(xi=0;xi<nx;xi++){
        for(yi=0;yi<ny;yi++){
            for(zi=0;zi<nz;zi++){
                density_field[zi+nz*(yi + ny*xi)] = 0. ;
                }
        }
     }
    /*calculate density field from particle info first*/
    for(i=1;i<=particle_number;i++)
	    {
            xi = (int)(P[i].Pos[0]/box_unit_length);
            yi = (int)(P[i].Pos[1]/box_unit_length);
            zi = (int)(P[i].Pos[2]/box_unit_length);

            l01 = P[i].Pos[0] - box_unit_length * (int)(P[i].Pos[0]/box_unit_length);
            l02 = box_unit_length - l01;
            l11 = P[i].Pos[1]- box_unit_length * (int)(P[i].Pos[1]/box_unit_length);
            l12 = box_unit_length - l11;
            l21 = P[i].Pos[2]- box_unit_length * (int)(P[i].Pos[2]/box_unit_length);
            l22 = box_unit_length - l21;
            
            max_weight = pow(pow((box_unit_length - maxoftwo(l01, l02)), 2) + pow((box_unit_length - maxoftwo(l11, l12)), 2) + pow((box_unit_length - maxoftwo(l21, l22)), 2), -1./2.);
            
            a_weight[0] = pow(l01*l01+l11*l11+l21*l21, -1./2.);
            a_weight[1] = pow(l02*l02+l11*l11+l21*l21, -1./2.);
            a_weight[2] = pow(l01*l01+l12*l12+l21*l21, -1./2.);
            a_weight[3] = pow(l02*l02+l12*l12+l21*l21, -1./2.);
            a_weight[4] = pow(l01*l01+l11*l11+l22*l22, -1./2.);
            a_weight[5] = pow(l02*l02+l11*l11+l22*l22, -1./2.);
            a_weight[6] = pow(l01*l01+l12*l12+l22*l22, -1./2.);
            a_weight[7] = pow(l02*l02+l12*l12+l22*l22, -1./2.);
            
            if(fractal_density == 0){
            
            if(a_weight[0]==max_weight){
                density_field[zi+nz*(yi + ny*xi)] += P[i].Mass;
            }
            
            if(a_weight[1]==max_weight){
                density_field[zi+nz*(yi + ny*(xi+1))] += P[i].Mass;
            }
            
            if(a_weight[2]==max_weight){
                density_field[zi+nz*((yi+1) + ny*xi)] += P[i].Mass;
            }
            
            if(a_weight[3]==max_weight){
                density_field[zi+nz*((yi+1) + ny*(xi+1))] += P[i].Mass;
            }
            
            if(a_weight[4]==max_weight){
                density_field[(zi+1)+nz*(yi + ny*xi)] += P[i].Mass;
            }
            
            if(a_weight[5]==max_weight){
                density_field[(zi+1)+nz*(yi + ny*(xi+1))] += P[i].Mass;
            }
            
            if(a_weight[6]==max_weight){
                density_field[(zi+1)+nz*((yi+1) + ny*xi)] += P[i].Mass;
            }
            
            if(a_weight[7]==max_weight){
                density_field[(zi+1)+nz*((yi+1) + ny*(xi+1))] += P[i].Mass;
            }
            }
            
            if(fractal_density == 1){
            //unit_weight = P[i].Mass * 1./(a_weight[0]+a_weight[1]+a_weight[2]+a_weight[3]+a_weight[4]+a_weight[5]+a_weight[6]+a_weight[7]);
            unit_weight = 1./(a_weight[0]+a_weight[1]+a_weight[2]+a_weight[3]+a_weight[4]+a_weight[5]+a_weight[6]+a_weight[7]);
            density_field[zi+nz*(yi + ny*xi)] += a_weight[0]*unit_weight;
            density_field[zi+nz*(yi + ny*(xi+1))] += a_weight[1]*unit_weight;
            density_field[zi+nz*((yi+1) + ny*xi)] += a_weight[2]*unit_weight;
            density_field[zi+nz*((yi+1) + ny*(xi+1))] += a_weight[3]*unit_weight;
            density_field[(zi+1)+nz*(yi + ny*xi)] += a_weight[4]*unit_weight;
            density_field[(zi+1)+nz*(yi + ny*(xi+1))] += a_weight[5]*unit_weight;
            density_field[(zi+1)+nz*((yi+1) + ny*xi)] += a_weight[6]*unit_weight;
            density_field[(zi+1)+nz*((yi+1) + ny*(xi+1))] += a_weight[7]*unit_weight;
            }
	 }
   
   printf("checkpoint4\n");
   
//----------------------------------------below is for distribution evolution correction called frstr--------------------

if(correction_nu_distribution == 1){

   FILE *fload_e;
   fload_e = fopen(f_load_early, "r");

   printf("checkpoint 4.1 %d\n", nnn);
  
  //load the delta k0 from txt 
   for(i=0;i<nnn;i++){
      fscanf(fload_e, "%lf ", &load_early_re[i]);
      //printf("%lf", load_early_re[i]);
   }

   

   FILE *fload_e2;
   fload_e2 = fopen(f_load_early2, "r");

      for(i=0;i<nnn;i++){
      fscanf(fload_e2, "%lf ", &load_early_im[i]);
      //printf("%lf ", load_early_re[i]);
   }

   for(i=0;i<nnn;i++){
    load_early[i][0] = load_early_re[i];
    load_early[i][1] = load_early_im[i];
    delta_ka1[i][1] = 0.;
    delta_ka1[i][0] = 0.;
    //printf("%lf %d\n", load_early[i][0], i);
   }

   printf("checkpoint 4.2\n");
 

   printf("checkpoint 4.5\n");
   //do the evolution for neutrino distribution in k space
   double kkk;

      //this part is the correction of neutrino spatial distribution
      printf("initialze array finished\n");
      
      
      for(xi=0;xi<=one_d_box;xi++){
            for(yi=0;yi<=one_d_box;yi++){
                //#pragma omp parallel for num_threads(4)
                for(zi=0;zi<=one_d_box;zi++){
                  //printf("zi %d\n", zi);
                  
                  if(xi>nmid){
                    xtemp = one_d_box + 1 - xi;
                  }
                  else{
                    xtemp =  - xi;
                  }

                  if(yi>nmid){
                    ytemp = one_d_box + 1 - yi;
                  }
                  else{
                    ytemp =  - yi;
                  }

                  if(zi>nmid){
                    ztemp = one_d_box + 1 - zi;
                  }
                  else{
                    ztemp =  - zi;
                  }
                 // #pragma omp critical
                  //{
                  kx = (2*M_PI*1000/box_s) * xtemp;  //MPc ^ -1
                  ky = (2*M_PI*1000/box_s) * ytemp;
                  kz = (2*M_PI*1000/box_s) * ztemp;
                
                  kkk = sqrt(kx*kx + ky*ky + kz*kz);
                //}
                  if(kkk > kcritical){
                    continue;
                  }
                  else{
                  //printf("how about this %lf %lf %lf %lf %lf\n", kkk, array_s[frstr_interval], array_s[0], kkk*array_s[frstr_interval], kkk*array_s[0]);
                  //printf("load early re %lf, load early im %lf\n", load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]);                 //
                 // #pragma omp critical
                  //  {
                  delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = frstr(kkk, load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], array_a, array_s);
                  delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] / load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]) * delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1];
                  //delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = frstr(kkk, load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], array_a, array_s);
                 // }
                  double r1 = delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] / load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0];
                  double r2 = delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] / load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1];
                  
                  //printf("re check for the integration solver %d %d %d delta k %.20f, k0 %lf\n", xi, yi, zi, delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]);
                  //printf("im check for the integration solver %d %d %d delta k %.20f, k0 %lf\n", xi, yi, zi, delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], load_early[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]);

                  if(r1<0.01 & r2 < 0.01){
                    break;
                  }
                }
                }
              }
            }
          printf("checkpoint 5\n");
    


    delta_ka1[0][0] = 0.;
    delta_ka1[0][1] = 0.;
   //F transform the neutrino fluctuation in k space
    fftw_execute(plan_of_nu);
    /*
      for(xi=0;xi<=one_d_box;xi++){
            for(yi=0;yi<=one_d_box;yi++){
                for(zi=0;zi<=one_d_box;zi++){
                  printf("delta_ro re %lf im %lf %d %d %d\n", delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], xi, yi, zi);
               
                  }
                }
              }
              */

    for(xi=0;xi<nx;xi++){
        for(yi=0;yi<ny;yi++){
            for(zi=0;zi<nz;zi++){
              delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] / pow((one_d_box+1), 3);
              //printf("delta ro xi yi zi %d %d %d, re %lf im %lf\n", xi, yi, zi, delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]);
            }
          }
        }

      }  

      //----------------------------------------above is for distribution evolution correction called frstr--------------------

    average = (double)(particle_number) / (double)(pow(one_d_box+1, 3));
    printf("mass = %f average value = %f\n", P[1].Mass, average);
    double a_scale = 1.;

    for(i=0;i<nx*ny*nz;i++){
        fluc_density_field[i][0] = 0.;
        fluc_density_field[i][1] = 0.;
        out[i][0] = 0.;
        out[i][1] = 0.;
    }
      //----------------------------------------below is direct perterbation with gravitational potential
    
      if(correction_nu_distribution == 2){
    for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
                //fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)] - average)/average;
              fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)])/average;
                fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = 0.;
                if(fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] < 0.){
                
                printf("dens = %f fluc = %f ave = %f xyz = %d %d %d\n", density_field[zi+nz*(yi + ny*xi)], fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)], average, xi, yi, zi);
                }
             }
        }
    }
    
    fftw_execute(fluc_density_field_plan);

/*
    fftw_complex *out_inver;

    out_inver = (fftw_complex*) malloc((one_d_box + 1)*(one_d_box + 1)*(one_d_box + 1) * sizeof(fftw_complex));

    fftw_plan test_out_inver;
    test_out_inver = fftw_plan_dft_3d((one_d_box+1), (one_d_box+1), (one_d_box+1), out, out_inver, -1, FFTW_MEASURE);

    fftw_execute(test_out_inver);


    for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
         
              printf("test imaginary %.20f real %.20f\n", out_inver[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][1], out_inver[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0]);
            }
          }
        }
*/

  fftw_complex *phik;
  fftw_complex *phi_real;

  phik = (fftw_complex*) malloc((one_d_box + 1)*(one_d_box + 1)*(one_d_box + 1) * sizeof(fftw_complex));
  phi_real = (fftw_complex*) malloc((one_d_box + 1)*(one_d_box + 1)*(one_d_box + 1) * sizeof(fftw_complex));

  fftw_plan plan_of_phi;
  plan_of_phi = fftw_plan_dft_3d((one_d_box+1), (one_d_box+1), (one_d_box+1), phik, phi_real, -1, FFTW_MEASURE);

  for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
         
              phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = 0.;
              phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = 0.;
              phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = 0.;
              phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = 0.;
                //printf("out %f %f\n", out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]);
            }
          }
        }


 
     for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){

                  if(xi>nmid){
                    xtemp = one_d_box + 1 - xi;
                  }
                  else{
                    xtemp =  - xi;
                  }

                  if(yi>nmid){
                    ytemp = one_d_box + 1 - yi;
                  }
                  else{
                    ytemp =  - yi;
                  }

                  if(zi>nmid){
                    ztemp = one_d_box + 1 - zi;
                  }
                  else{
                    ztemp =  - zi;
                  }

                  kx = (2*M_PI*1000/box_s) * xtemp;  //MPc ^ -1
                  ky = (2*M_PI*1000/box_s) * ytemp;
                  kz = (2*M_PI*1000/box_s) * ztemp;
                
                  k = sqrt(kx*kx + ky*ky + kz*kz);

                  if(k == 0){
                    phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = 0.;
                    phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = 0.;
                  }

                  else{
                  phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (0. - out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]) * rocr / cc * Omega_m * 4. * M_PI * Gr / ((a_scale * a_scale * a_scale) * (k * k / mpc_to_m / mpc_to_m) * pow((one_d_box+1), 3));
                  phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = (0. - out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1]) * rocr / cc * Omega_m * 4. * M_PI * Gr / ((a_scale * a_scale * a_scale) * (k * k / mpc_to_m / mpc_to_m) * pow((one_d_box+1), 3));
                  //printf("phi = %lf, k = %lf out = %lf norm %f\n", phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], mpc_to_m * mpc_to_m / k / k / 1e40, out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] * rocr * Omega_m * 4. * M_PI * Gr * 1e15, pow((one_d_box+1), 3));
                }
            }
          }
        }

        
        fftw_execute(plan_of_phi);


     for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){

           //phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0];
           //phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1];
           //phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] = phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1] / pow((one_d_box+1), 3); 
           //printf("phi_real imaginary %f real %f phik %lf\n", phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], phik[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]);
           //printf("phi/T %f phi %f\n", phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] * 0.1 / ktoev / cc / 1.7, phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]);
            }
          }
        }

        double ev_to_kg = 1.78266191E-36;
        double mass_of_sun = 1.989E30;
        double average2;
        double rosum = 0.;
        average2 = rocr * Omega_m * pow((200 * mpc_to_m / h), 3) / cc / mass_of_sun / 1e10 * h / 33.5;

      for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
           
              delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = mass_neu * ev_to_kg * perturb_fd_integration(a_scale, phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]) * pow((mpc_to_m * box_s / 1000. / (one_d_box+1)), 3) / h / h / mass_of_sun / 1E10 / P[1].Mass;
              
              rosum = rosum + delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] / average;
              //delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = mass_neu * ev_to_kg * perturb_fd_integration(a_scale, 0.) * pow((mpc_to_m * box_s / 1000. / one_d_box), 3) * h * h / mass_of_sun / 1E10 / P[1].Mass;
              //delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = perturb_fd_integration(a_scale, 0.) * pow((mpc_to_m * box_s / 1000. / one_d_box), 3) * h * h / mass_of_sun / 1E10 / P[1].Mass / cc;
              //double cdm = density_field[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)];
              //double neuneu = mass_neu * ev_to_kg * perturb_fd_integration(a_scale, phi_real[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]) * cc;
              
              /*if(cdm > 1){
                 printf("delta_ro = %lf, cdm mass density %lf mass %lf average %lf ro %lf omega ratio %lf\n", delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0] / average, cdm, P[1].Mass, average, Omega_nu0);
              }*/
             
            }
          }
        }
        rosum = rosum / pow(one_d_box+1, 3);
        printf("rosum/average %lf Omega_nu0 / Omega_m %lf\n", rosum, neutrino_integration(1.) / rocr / Omega_m);

    }

    //-------------------------------------potential correction finished-----------------------------

    for(xi=0;xi<4;xi++){
        for(yi=0;yi<4;yi++){
            for(zi=0;zi<4;zi++){
              printf("delta ro xi yi zi %d %d %d, re %lf im %lf delta_ka1 %f %f\n", xi, yi, zi, delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], delta_ka1[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0]);
            }
          }
        }

    double delta_ro_sum = 0.;

    for(xi=0;xi<nx;xi++){
        for(yi=0;yi<ny;yi++){
            for(zi=0;zi<nz;zi++){
              delta_ro_sum = delta_ro_sum + delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0];
              //printf("delta_ro_sum %lf\n", delta_ro_sum);
            }
          }
        }

        printf("--------------sum of delta ro = %f ----------------\n", delta_ro_sum);
    
    /*delta density field and fourier transform it*/

    
    for(i=0;i<nx*ny*nz;i++){
        fluc_density_field[i][0] = 0.;
        fluc_density_field[i][1] = 0.;
        out[i][0] = 0.;
        out[i][1] = 0.;
    }
    
     /*calculate the fluctuation fueld and execute the fourier transform plan*/
    printf("look here\n");
    for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
                //fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)] = (density_field[zi+nz*(yi + ny*xi)] - average)/average;
                //printf("look here again\n");
                double roro = delta_ro[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0];
                if(correction_nu_distribution == 1){
                //fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)] * Omega_m / (Omega_m + Omega_nu0) + (roro + 1.) * average * Omega_nu0 / (Omega_m + Omega_nu0) - average)/average;
                fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)] + (roro + 1.) * average * Omega_nu0 / (Omega_m) - average)/average;
                }

                if(correction_nu_distribution == 2){
                fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)] * Omega_m / (Omega_m + Omega_nu0) + roro - average)/average; 
                }
                //if(fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)] != 0.){
                
                //printf("dens = %f fluc = %f ro %lf ave = %f xyz = %d %d %d\n", density_field[zi+nz*(yi + ny*xi)], fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0], roro, average, xi, yi, zi);
                //}
             }
        }
    }

    //printf("check a\n");

    fftw_execute(fluc_density_field_plan);

    /*    
      for(xi=0;xi<=one_d_box;xi++){
            for(yi=0;yi<=one_d_box;yi++){
                for(zi=0;zi<=one_d_box;zi++){
                  printf("out re %lf im %lf %d %d %d\n", out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][0], out[zi+(one_d_box+1)*(yi + (one_d_box+1)*xi)][1], xi, yi, zi);
               
                  }
                }
              }*/
              
    fftw_free(density_field);
    
    /*count the number of k and pk points in each bin*/
   for(b=0;b<Nbins;b++){
       for(xi=0;xi<=one_d_box;xi++){
        for(yi=0;yi<=one_d_box;yi++){
            for(zi=0;zi<=one_d_box;zi++){
                  if(xi>nmid){
                    xtemp = one_d_box + 1 - xi;
                  }
                  else{
                    xtemp =  - xi;
                  }

                  if(yi>nmid){
                    ytemp = one_d_box + 1 - yi;
                  }
                  else{
                    ytemp =  - yi;
                  }

                  if(zi>nmid){
                    ztemp = one_d_box + 1 - zi;
                  }
                  else{
                    ztemp =  - zi;
                  }

                  kx = (2*M_PI*1000/box_s) * xtemp;  //MPc ^ -1
                  ky = (2*M_PI*1000/box_s) * ytemp;
                  kz = (2*M_PI*1000/box_s) * ztemp;
                
                  k = sqrt(kx*kx + ky*ky + kz*kz);
                    if(k>kbins0[b] && k<=kbins0[b+1]){
                        count[b]++;
                    }
                }
            }
        }
   //printf("kb = %f number of kb = %d\n", kbins[b], count[b]);
   }
   //printf("check b\n");
    //print to check 3d fourier box
    /*
    double pkkk;
    for(xi=0;xi<nx;xi++){
        for(yi=0;yi<ny;yi++){
            for(zi=0;zi<nz;zi++){
                kx = (2*M_PI*1000/box_s) * xi;
                ky = (2*M_PI*1000/box_s) * yi;
                kz = (2*M_PI*1000/box_s) * zi;
                
                double rr = out[zi+nz*(yi + ny*xi)][0];
                double ii = out[zi+nz*(yi + ny*xi)][1];
                k = sqrt(kx*kx + ky*ky + kz*kz);
                    pkkk = rr*rr + ii*ii;
                if(rr ==0. && ii != 0.){
                    printf("here real part is 0 %d %d %d\n", xi, yi, zi);
                }
                if(ii ==0 && rr != 0.){
                    printf("here img part is 0 %d %d %d\n", xi, yi, zi);
                }
               // printf("xyz = %d %d %d    rr = %f   ii = %f  pk = %f\n", xi, yi, zi, rr, ii, pkkk);
                if(ii!=0 || rr!=0){
                
                }
                }
            }
        }
    printf("00r = %f  00i = %f\n", out[zi+nz*(yi + ny*xi)][0], out[zi+nz*(yi + ny*xi)][1]);
    */
    /*use the counted number of k points to calculate sum, average and std in each bin*/
    for(b=0;b<Nbins;b++){
        double *kdata;
        kdata = (double*) malloc(count[b] * sizeof(double));
        double pkk;
        double sum, var, std, ave;
        kerr[b] = 0.;
        pk[b] = 0.;
        
        for(i=0;i<count[b];i++){
            kdata[i] = 0.;
        }
        
        for(xi=0;xi<=one_d_box;xi++){
            for(yi=0;yi<ny;yi++){
                for(zi=0;zi<nz;zi++){
                  if(xi>nmid){
                    xtemp = one_d_box + 1 - xi;
                  }
                  else{
                    xtemp =  - xi;
                  }

                  if(yi>nmid){
                    ytemp = one_d_box + 1 - yi;
                  }
                  else{
                    ytemp =  - yi;
                  }

                  if(zi>nmid){
                    ztemp = one_d_box + 1 - zi;
                  }
                  else{
                    ztemp =  - zi;
                  }

                  kx = (2*M_PI*1000/box_s) * xtemp;  //MPc ^ -1
                  ky = (2*M_PI*1000/box_s) * ytemp;
                  kz = (2*M_PI*1000/box_s) * ztemp;                  

                   k = sqrt(kx*kx + ky*ky + kz*kz);

                    if(k>kbins0[b] && k<=kbins0[b+1]){
                       // int nu;
                       // pkk = (out[zi+nz*(yi + ny*xi)][1])*(out[zi+nz*(yi + ny*xi)][1]) + (out[zi+nz*(yi + ny*xi)][0])*(out[zi+nz*(yi + ny*xi)][0]);
                         pkk = (out[zi+nz*(yi + ny*xi)][1])*(out[zi+nz*(yi + ny*xi)][1]) + (out[zi+nz*(yi + ny*xi)][0])*(out[zi+nz*(yi + ny*xi)][0]);
                        //printf("k = %f    powerk = %f xi %d yi %d zi %d\n", k, pkk, xi, yi, zi);
                        
                        kdata[count1[b]] = pkk;
                        count1[b]++;
                    }
                }
            }
        }
        //printf("kb = %f number of kb = %d\n", kbins[b], count1[b]);
        
        if(count[b]!= 0){
            sum = 0.0;
            for(i=0;i<count[b];i++){
                sum += kdata[i];
            }
            var = 0.0;
            std = 0.;
            ave = sum/count[b];
            pk[b] = ave*pk_nor;
            //printf("pk = %f\n", pk[b]);
            for(i=0;i<count[b];i++){
                var += (kdata[i] - ave)*(kdata[i] - ave);
            }
            std = sqrt(var);
            kerr[b] = std*pk_nor/count[b];
            //printf("kbin = %f count[b] = %d pk[b] = %lf ave = %f err = %f\n", kbins[b], count[b], pk[b], ave, kerr[b]);
        }

        pk_corrected[b] = pk[b];
        kerr_corrected[b] = kerr[b];

        //doing the yp jing correction of iteration

        double prk, alpha0, alpha1, alpha2, sum2, sum3, yuan, tta, c2, p0;
        iteration_count = 0;
        if(count[b]!= 0){
        alpha0 = log(pk[b]) / log(kbins[b]/10.);
        alpha1 = alpha0 + 1.;
        alpha2 = alpha0;
        prk = pk[b] - (1. - (2./3.) * pow((sin(M_PI * kbins[b] / (2 * M_PI * 1000./box_unit_length))), 2))/ particle_number;

        //double testnum = w(2*kn + (2*M_PI/box_s), 2*kn + (2*M_PI/box_s), 2*kn + (2*M_PI/box_s));
        //printf("?? %f testnum%15f\n", (2*M_PI/box_s), w(2*kn + (2*M_PI/box_s), 2*kn + (2*M_PI/box_s), 2*kn + (2*M_PI/box_s)));
        //printf("iteration %d alpha0 %f, alpha1 %f, alpha2 %f\n", iteration_count, alpha0, alpha1, alpha2);
        //if(kbins[b] > kn*1000./8.){
        if(kbins[b] > 0.2){
        while (fabs(alpha2 - alpha1) > (fabs(alpha1) * 0.001)){
          sum3 = 0.;
          alpha1 = alpha2;
          for(xi=1;xi<=one_d_box;xi++){
            for(yi=1;yi<ny;yi++){
                for(zi=1;zi<nz;zi++){
                  if(xi>nmid){
                    xtemp = one_d_box + 1 - xi;
                  }
                  else{
                    xtemp =  - xi;
                  }

                  if(yi>nmid){
                    ytemp = one_d_box + 1 - yi;
                  }
                  else{
                    ytemp =  - yi;
                  }

                  if(zi>nmid){
                    ztemp = one_d_box + 1 - zi;
                  }
                  else{
                    ztemp =  - zi;
                  }

                  kx = (2*M_PI*1000/box_s) * xtemp;  //MPc ^ -1
                  ky = (2*M_PI*1000/box_s) * ytemp;
                  kz = (2*M_PI*1000/box_s) * ztemp;
                  
                  k = sqrt(kx*kx + ky*ky + kz*kz);

                    if(k>kbins0[b] && k<=kbins0[b+1]){
                      sum2 = 0.;
                      for(i=0;i<5;i++){
                        for(j=0;j<5;j++){
                          for(k=0;k<5;k++){
                            //tta = pow((kx + 2*kn*1000.*i), 2) + pow((ky + 2*kn*1000.*j), 2) + pow((kz + 2*kn*1000.*k), 2);
                            tta = pow((kx/10. + 2*kn*100*i), 2) + pow((ky/10. + 2*kn*100*j), 2) + pow((kz/10. + 2*kn*100*k), 2);
                        yuan = pow(w((kx/1000. + 2*kn*i), (ky/1000. + 2*kn*j), (kz/1000. + 2*kn*k)), 2) * pow(tta, alpha1/2.);
                        sum2 = sum2 + yuan;
                        //if (pow(tta, alpha1/2.) > 0.001){
                        //printf("kxx %f kyy %f kzz %f inter %f \n", (kx + 2*kn*1000.*i), (ky + 2*kn*1000.*j), (kz + 2*kn*1000.*k), M_PI*(kz/1000. + 2*kn*k)/(2.*kn));
                        //printf("k %f pkb %f, k1 %f, tta %f power %f, powerpk %f, yuan %f\n", kbins[b], pk[b], (kx/1000. + 2*kn*i), tta, alpha1/2., pow(tta, alpha1/2.), yuan);
                        
                      }
                    //}
                }
               }
           sum3 = sum3 + sum2; }
          }
        
        } 
      }
      c2 = (sum3/count[b])/pow(kbins[b]/10., alpha1);
      p0 = prk / c2;
      alpha2 = log(p0) / log(kbins[b]/10.);
      iteration_count = iteration_count + 1;
      //printf("measure %f , calculate %f k %f\n", pk[b], sum3/count[b], kbins[b]);
      //printf("iteration %d alpha0 %f, alpha1 %f, alpha2 %f delta %f \n", iteration_count, alpha0, alpha1, alpha2, fabs(alpha2 - alpha1));
      //sleep(1);
    }
    pk_corrected[b] = pow(kbins[b]/10., alpha2);
  }
          
        }

    }
    //fftw_destroy_plan(fluc_density_field_plan);

    /*for(b=0;b<Nbins;b++){
      printf("kbin %lf, pk_un %lf, pk %lf\n", kbins[b], pk[b], pk_corrected[b]);
    }*/

    char output[200];
    //sprintf(output, "%s/frstr%.2f.txt", addfload, frstr_a0);
    sprintf(output, "%s/%s.txt", addfload, txtname);
    FILE *fp;
    fp=fopen(output,"w");
    for(b=0;b<Nbins;b++) {
        fprintf(fp,"%f\t %f\t %f\t %f\n", kbins[b], pk_corrected[b], pk[b], kerr[b]);
    }
    fclose(fp);
    
    fftw_free(out);
    fftw_free(load_early);
    fftw_free(delta_ka1);
    fftw_free(fluc_density_field);
    fftw_free(delta_ro);
    
    return 0;
}





/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(void)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;

  /* physical constants in cgs units */
  GRAVITY = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm = 3.085678e21;	/*  code length unit in cm/h */
  UnitMass_in_g = 1.989e43;	/*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s = 1.0e5;

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);


  Xh = 0.76;			/* mass fraction of hydrogen */
  HubbleParam = 0.65;


  for(i = 1; i <= NumPart; i++)
    {
      if(P[i].Type == 0)	/* gas particle */
	{
	  MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * P[i].Ne) * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u = P[i].U * UnitEnergy_in_cgs / UnitMass_in_g;

	  gamma = 5.0 / 3;

	  /* get temperature in Kelvin */

	  P[i].Temp = MeanWeight / BOLTZMANN * (gamma - 1) * u;
	}
    }
    printf("unittrans is ok\n");
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, j, k, dummy, ntot_withmasses;
  int t, n, off, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading `%s' ...\n", buf);
      fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npart[k];
	  Ngas = header1.npart[0];
	}
      else
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npartTotal[k];
	  Ngas = header1.npartTotal[0];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	allocate_memory();

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	       pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      P[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1.npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  Time = header1.time;
  Redshift = header1.time;
  printf("pc_new = %d\n", pc_new);
  particle_number = pc_new - 1;
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P = malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  P--;				/* start with offset 1 */


  if(!(Id = malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  Id--;				/* start with offset 1 */

  printf("allocating memory...done\n");
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i, j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i = 1; i <= NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource = P[i];
	  idsource = Id[i];
	  dest = Id[i];

	  do
	    {
	      psave = P[dest];
	      idsave = Id[dest];

	      P[dest] = psource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      psource = psave;
	      idsource = idsave;

	      dest = idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;
  free(Id);

  printf("space for particle ID freed\n");
}


