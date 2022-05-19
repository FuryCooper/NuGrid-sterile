#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include "vars.h"

//this version the unit is 1/Mpc instead of h/Mpc
// and measures \Delta^2 = Pk * pow(k, 3)

int particle_number, one_d_num, one_d_box;
double box_s;
double hubble;
#define M_PI 3.14159265358979323846
#define Nbins 28
#define field_resolution 128
#define kbininterval 1.14
#define start_k 0.04

int shift_on = 0;

int fractal_density = 1;
int reading_today = 1;
double pk_nor;      //normalization constant for power spectrum

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

    char path[300], input_fname[300], basename[300];
  int type, snapshot_number, files, kk;
  int start_time, finish_time;
    
    read_parameterfile(argv[1]);
    
    snapshot_number = snapno;		/* number of snapshot */
    files = 1;			/* number of files per snapshot */

    if (reading_today == 1){
        sprintf(path, input_path);
         sprintf(basename, "snapshot");
        sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
    }
    else {sprintf(path, "/home/zzc/ICs/mar1");
        sprintf(basename, "123457-random-ics");
        sprintf(input_fname, "%s/%s", path, basename);
    }
    
  start_time = clock();
  load_snapshot(input_fname, files);

   box_s = header1.BoxSize;
    hubble = header1.HubbleParam;
    
    //Omega_lambda = header1.OmegaLambda;
    //Omega_m = header1.Omega0;
    
  printf("particle number is %d mass[1] = %f hubbleparam = %f\n", particle_number, P[1].Mass, hubble);
  pk_nor = pow(box_s/1e3, 3)/(pow(field_resolution+1, 6));
  one_d_num = (int)(pow(particle_number, 1./3.)+0.5);

  one_d_box = field_resolution;
  printf("Separating the simulation box into %d grids...\n", one_d_box);
    
  reordering();			/* call this routine only if your ID's are set properly */
  unit_conversion();		/* optional stuff */
  //printf("Mass3 = %f\n", P[3].Mass);
  do_what_you_want();
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
    double up1, up2, up3, down;
    double kn = M_PI / (box_s / (double)(field_resolution));
    
    if(fabs(k1)<1e-3){
        up1 = 1. - pow(M_PI * k1 / (2. * kn), 2) / 6.;
    }
    
    else{
        up1 = sin(M_PI * k1 / (2. * kn)) / (M_PI * k1 / (2. * kn));
    }
    
    if(fabs(k2)<1e-3){
        up2 = 1. - pow(M_PI * k2 / (2. * kn), 2) / 6.;
    }
    
    else{
        up2 = sin(M_PI * k2 / (2. * kn)) / (M_PI * k2 / (2. * kn));
    }
    
    if(fabs(k3)<1e-3){
        up3 = 1. - pow(M_PI * k3 / (2. * kn), 2) / 6.;
    }
    
    else{
        up3 = sin(M_PI * k3 / (2. * kn)) / (M_PI * k3 / (2. * kn));
    }
    
    down = (up1*up2*up3)*(up1*up2*up3);
    return down;
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
    double pk[Nbins], pk_corrected[Nbins], kerr_corrected[Nbins], pk_win[Nbins];
    
    double kn = M_PI / (box_s / (double)(field_resolution));
    
    fluc_density_field = (fftw_complex*) fftw_malloc((one_d_box+1)*(one_d_box+1)*(one_d_box+1) * sizeof(fftw_complex));
    out = (fftw_complex*) fftw_malloc(nx*ny*nz * sizeof(fftw_complex));

    fftw_plan fluc_density_field_plan;
    fluc_density_field_plan = fftw_plan_dft_3d((one_d_box+1), (one_d_box+1), (one_d_box+1), fluc_density_field, out, 1, FFTW_MEASURE);

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
        pk_corrected[i] = 0.;
        pk_win[i] = 0.;
    }
    
    for(i=0;i<Nbins;i++){
        //kbins[i] = (kbins0[i]+kbins1[i])/2.;
        kbins[i] = sqrt(kbins0[i]*kbins1[i]);
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

    /*calculate density field from particle info mesh*/
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
        /*
        max_weight = pow(pow((box_unit_length - maxoftwo(l01, l02)), 2) + pow((box_unit_length - maxoftwo(l11, l12)), 2) + pow((box_unit_length - maxoftwo(l21, l22)), 2), -1./2.);
         
         a_weight[0] = pow(l01*l01+l11*l11+l21*l21, -1./2.);
         a_weight[1] = pow(l02*l02+l11*l11+l21*l21, -1./2.);
         a_weight[2] = pow(l01*l01+l12*l12+l21*l21, -1./2.);
         a_weight[3] = pow(l02*l02+l12*l12+l21*l21, -1./2.);
         a_weight[4] = pow(l01*l01+l11*l11+l22*l22, -1./2.);
         a_weight[5] = pow(l02*l02+l11*l11+l22*l22, -1./2.);
         a_weight[6] = pow(l01*l01+l12*l12+l22*l22, -1./2.);
         a_weight[7] = pow(l02*l02+l12*l12+l22*l22, -1./2.);
        
        */
        max_weight = maxoftwo(l01, l02) * maxoftwo(l11, l12) * maxoftwo(l21, l22) / pow(box_unit_length, 3);
        
        a_weight[0] = (1.-l01/box_unit_length)*(1.-l11/box_unit_length)*(1.-l21/box_unit_length);
        a_weight[1] = (1.-l02/box_unit_length)*(1.-l11/box_unit_length)*(1.-l21/box_unit_length);
        a_weight[2] = (1.-l01/box_unit_length)*(1.-l12/box_unit_length)*(1.-l21/box_unit_length);
        a_weight[3] = (1.-l02/box_unit_length)*(1.-l12/box_unit_length)*(1.-l21/box_unit_length);
        a_weight[4] = (1.-l01/box_unit_length)*(1.-l11/box_unit_length)*(1.-l22/box_unit_length);
        a_weight[5] = (1.-l02/box_unit_length)*(1.-l11/box_unit_length)*(1.-l22/box_unit_length);
        a_weight[6] = (1.-l01/box_unit_length)*(1.-l12/box_unit_length)*(1.-l22/box_unit_length);
        a_weight[7] = (1.-l02/box_unit_length)*(1.-l12/box_unit_length)*(1.-l22/box_unit_length);
        
        
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
            //unit_weight = 1. / (a_weight[0]+a_weight[1]+a_weight[2]+a_weight[3]+a_weight[4]+a_weight[5]+a_weight[6]+a_weight[7]);
            unit_weight = P[i].Mass;
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

    int xtemp, ytemp, ztemp;
    int nmid = one_d_box/2;
    
    average = 0.;
    
    for(i=0;i<particle_number;i++){
        average += P[i].Mass;
    }
    
    average *= 1. / (double)(pow(one_d_box+1, 3));
    
    printf("mass = %f average value = %f\n", P[1].Mass, average);
    printf("mass2 = %f average value = %f\n", P[particle_number-1].Mass, average);
    
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
                fluc_density_field[zi + (one_d_box+1)*(yi + (one_d_box+1)*xi)][0] = (density_field[zi+nz*(yi + ny*xi)] - average)/average;
                  }
        }
    }
    
    fftw_execute(fluc_density_field_plan);
    
    
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
                    if(k>kbins0[b] && k<=kbins1[b]){
                        count[b]++;
                    }
                }
            }
        }
   //printf("kb = %f number of kb = %d\n", kbins[b], count[b]);
   }
   
    
    //Now read the neutrino txt
    
    
    int ii;
    
    int rd, rd_size;
    
    FILE *fd;
    double kkk, ppp;
    double *rd_array_k, *rd_array_pk;
    
    rd_size = 0;
    
    if(!(fd = fopen(nu_txt, "r")))
    {
        printf("can't read input spectrum in file '%s'\n", nu_txt);
    }
    
    do
    {
        if(fscanf(fd, " %lg %lg ", &kkk, &ppp) == 2){
            rd_size++;
        }
        else
            break;
    }
    while(1);
    fclose(fd);
    
    rd_array_k = (double*) malloc((rd_size) * sizeof(double));
    rd_array_pk = (double*) malloc((rd_size) * sizeof(double));
    
    rd_size = 0;
    
    if(!(fd = fopen(nu_txt, "r")))
    {
        printf("can't read input spectrum in file '%s'\n", nu_txt);
    }
    
    do
    {
        if(fscanf(fd, " %lg %lg ", &kkk, &ppp) == 2)
        {
            rd_array_k[rd_size] = kkk;
            rd_array_pk[rd_size] = ppp;
            rd_size++;
        }
        else
            break;
    }
    while(1);
    fclose(fd);

    
    /*use the counted number of k points to calculate sum, average and std in each bin*/
    for(b=0;b<Nbins;b++){
        double *kdata, *kdata_win;
        kdata = (double*) malloc(count[b] * sizeof(double));
        kdata_win = (double*) malloc(count[b] * sizeof(double));
        double pkk, pkwin;
        double sum, var, std, ave, sum_win;
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

                    if(k>kbins0[b] && k<=kbins1[b]){
                         pkk = (out[zi+nz*(yi + ny*xi)][1])*(out[zi+nz*(yi + ny*xi)][1]) + (out[zi+nz*(yi + ny*xi)][0])*(out[zi+nz*(yi + ny*xi)][0]);
                        kdata[count1[b]] = pkk;
                        count1[b]++;
                    }
                }
            }
        }
        //printf("kb = %f number of kb = %d\n", kbins[b], count1[b]);
        
        
        
        if(count[b]!= 0){
            sum = 0.0;
            sum_win = 0.0;
            for(i=0;i<count[b];i++){
                sum += kdata[i];
            }
            var = 0.0;
            std = 0.;
            ave = sum/count[b];
            pk[b] = ave*pk_nor;
            //pk_win[b] = sum_win*pk_nor/count[b];
            //printf("pk = %f k = %f\n", pk[b], kbins[b]);
            for(i=0;i<count[b];i++){
                var += (kdata[i] - ave)*(kdata[i] - ave);
                //printf("kdata[i] %lf, ave %lf, ratio %lf\n", kdata[i], ave, kdata[i] / ave);
            }
            std = sqrt(var);
            kerr[b] = (std*pk_nor/ count[b]);
            //printf("kbin = %f count[b] = %d pk[b] = %lf ave = %f err = %f hubble %lf\n", kbins[b], count[b], pk[b], ave, kerr[b], hubble);
        }
        
        /*FILE *ddis;
        char dis_fname[200];
        int dd;
        sprintf(dis_fname, "Delta-dis-bin-%d.txt", b);
        ddis=fopen(dis_fname,"w");
        for(dd=0;dd<count[b];dd++) {
            fprintf(ddis,"%f\n", kdata[dd]);
        }
        fclose(ddis);*/

        
        double temp_ratio0, temp_ratio1, temp_ratio;
        
        if(count[b]!= 0){
            for(ii=0;ii<rd_size;ii++){
                
                if(kbins[b] >= (rd_array_k[ii]) && kbins[b] < (rd_array_k[ii+1])){
                    
                        temp_ratio0 = rd_array_pk[ii];
                    
                        temp_ratio1 = rd_array_pk[ii+1];
                    
                    temp_ratio = (temp_ratio0 + temp_ratio1) / 2.;
                }

            }
            
            if(kbins[b] > rd_array_k[rd_size]){
                temp_ratio = 0.;
            }
            
            pk[b] = pow(sqrt(pk[b]) * (1. - fnu) + sqrt(pk[b]) * sqrt(temp_ratio) * fnu, 2.);
            kerr[b] = pow(sqrt(kerr[b]) * (1. - fnu) + sqrt(kerr[b]) * sqrt(temp_ratio) * fnu, 2.);
        }
        
        
        pk_corrected[b] = pk[b];
        kerr_corrected[b] = kerr[b] / pk[b];

        //doing the yp jing correction of iteration

        double prk, alpha0, alpha1, alpha2, sum2, sum3, yuan, tta, c2, p0, clog;
        iteration_count = 0;
        if(b>0){
            if(count[b]!= 0 && count[b-1]!=0){
        alpha0 = (log(pk[b-1]) - log(pk[b]))/ (log(kbins[b-1]/10.) - log(kbins[b]/10.));
        clog = pk[b-1] / pow(kbins[b-1], alpha0);
        alpha1 = alpha0 + 1.;
        alpha2 = alpha0;
        prk = pk[b] - (1. - (2./3.) * pow((sin(M_PI * kbins[b] / (2 * M_PI * 1000./box_unit_length))), 2))/ particle_number;

        while (fabs(alpha2 - alpha1) > (fabs(alpha1) * 0.001)){
          sum3 = 0.;
          alpha1 = alpha2;
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

                    if(k>kbins0[b] && k<=kbins1[b]){
                      sum2 = 0.;
                      for(i=0;i<5;i++){
                        for(j=0;j<5;j++){
                          for(k=0;k<5;k++){
                            tta = pow((kx/10. + 2*kn*100*i), 2) + pow((ky/10. + 2*kn*100*j), 2) + pow((kz/10. + 2*kn*100*k), 2);
                        yuan = pow(w((kx/1000. + 2*kn*i), (ky/1000. + 2*kn*j), (kz/1000. + 2*kn*k)), 2) * pow(tta, alpha1/2.)*clog;
                        sum2 = sum2 + yuan;
                    }
                }
               }
           sum3 = sum3 + sum2; }
          }
        
        } 
      }
      c2 = (sum3/count[b])/pow(kbins[b]/10., alpha1) / clog;
      p0 = prk / c2;
      alpha2 = log(p0/clog) / log(kbins[b]/10.);
      iteration_count = iteration_count + 1;
    }
    pk_corrected[b] = clog*pow(kbins[b]/10., alpha2);
  }
          
        }

    }

    FILE *fp;
    fp=fopen(output,"w");
    for(b=0;b<Nbins;b++) {
        fprintf(fp,"%f\t %f\t %f\t %f %d\n", kbins[b], pk_corrected[b], pk[b], kerr_corrected[b], count[b]);
    }
    fclose(fp);
    
    fftw_free(out);
    fftw_free(fluc_density_field);
    fftw_free(density_field);
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
  char buf[300];
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
    
    for(n=1;n<20;n++){
        printf("id %d mass %f vel1 %f vel 2 %f vel3 %f\n", n, P[n].Mass, P[n].Vel[0], P[n].Vel[1], P[n].Vel[2]);
        
        printf("id %d mass %f vel1 %f vel 2 %f vel3 %f\n", n+128*128*128, P[n+128*128*128].Mass, P[n+128*128*128].Vel[0], P[n+128*128*128].Vel[1], P[n+128*128*128].Vel[2]);
    }
    
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
