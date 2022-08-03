/**************************************************************************
 
 Copyright (c) 2022 Neil Cornish
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "SpecUtils.h"
#include "Constants.h"

#ifndef _OPENMP
#define omp ignore
#endif


//##############################################
//OPEN MP

gsl_rng **rvec;
//##############################################


//OSX
// clang -O3 -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c SpecUtils.c -lgsl  -lm

// Linux
// gcc -O3 -std=gnu99 -fopenmp -w -o SpecFit SpecFit.c SpecUtils.c -lgsl -lm

void update(int m, int q, double *logLx, struct Smooth *smoothline, double heat, double smsc, double lnsc, double Tobs, int Ns, int Nknot, int Nlines, double *freqs, double *PS, double *LarrayX, double *SM, double *SL, double *SN, double *ffitx, double *Xsline,  double *bx, double *cx, double *slopex, double *delx,  double  *linef, double *lineh, double *lineQ, double *linew, double *deltafmax, int *cS, int *cL, int *acS, int *acL, gsl_rng *r);

int smoothset(int Ns, double *SM, double *freqs, double Tobs, double *fit, double *smline, int *Nk);

int lineset(int Ns, double *SM, double *PS, double *freqs, double Tobs, double *lf, double *lh, double *lQ, double *lw, double *dfmx, int *Nlns);

int main(int argc, char *argv[])
{

  int i, j, k, ii, kk, sc, ND, N, Ns, Nm, m, mc, dec, MM, MK, MS, skip;
  int scount, sacc, hold;
  int NC;
  double **LarrayX;
  double *logLx;
  int *who;
  double *heat;
  double **ffitx, **Xsline;
  int imin, imax;
  int *acS, *acL, *cS, *cL;
  int typ, sm, sma;
  double *timeF, *dataF;
  double *times, *data, *Hf;
  double x, y, z, dt, df, Tobs, ttrig, fny, fmn, fmx;
  double f, finc;
  double **SN, *SNA, *SMA, **SM, **SL, *PS, *sp;
    
    double lmax;
    int qmax;
 
  int Nsp, Nl;
  char command[1024];
  double spread, **deltafmax, **linew;
  double **linef, **lineh, **lineQ;
  double *smline, *fit;
  double *lf, *lh, *lQ, *lw, *dfmx;
  double f2, f4, ff4, ff2, df2, df4;
  int flag, Nlines;
  double xold, xnext, Abar;
  double max;
  int Nknot, Nsten;
  double a1, a2, sdd;
  double smsc, lnsc;
    
  // these are used by smooth linear fit. Have to be declared even if not allocated
  double **slopex, **ax, **cx, **bx, **delx;
  struct Smooth *smoothline  = malloc(sizeof(struct Smooth));
    
  // these are used by the spline. Have to be declared even if not allocated
  gsl_spline   *aspline;
  gsl_interp_accel *acc;
    
  double alpha, beta;
  double tuke = 0.4;    // Tukey window rise (s)
    
   const gsl_rng_type * T;
   gsl_rng * r;
    
   gsl_rng_env_setup();
    
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);
    
    clock_t start, end;
    double cpu_time_used;
    
    double itime, ftime, exec_time;
    
    
  FILE *in;
  FILE *out;
  FILE *spec;
  FILE *lfile;
  FILE *sfile;
    
   
    
    smsc = 1.0;
    lnsc = 1.0;
    
    if(argc!=5)
    {
        printf("./Spec Tobs trig_time obs fmax\n");
        return 1;
    }
    
    Tobs = atof(argv[1]);
    ttrig = atof(argv[2]);
    m = atoi(argv[3]);
    fmx = atof(argv[4]);
    
    df = 1.0/Tobs;
    
    MS = 1000;

    
    sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    
    in = fopen(command,"r");
 
    
    ND = -1;
    while(!feof(in))
    {
        fscanf(in,"%lf%lf", &x, &y);
        ND++;
    }
    rewind(in);
    printf("Number of points = %d\n", ND);
    
    timeF = (double*)malloc(sizeof(double)* (ND));
    dataF = (double*)malloc(sizeof(double)* (ND));
    
    for (i = 0; i < ND; ++i)
    {
      fscanf(in,"%lf%lf", &timeF[i], &dataF[i]);
    }
    fclose(in);
    
    dt = timeF[1]-timeF[0];
    Tobs = (double)(ND)*dt;  // duration
    fny = 1.0/(2.0*dt);  // Nyquist
    
    printf("Nyquist %f  fmax %f\n", fny, fmx);
    
    // if fmax < fny we can downsample the data by decimation
    // first we bandpass then decimate
    
    dec = (int)(fny/fmx);
    
    if(dec > 8) dec = 8;
    
    
    printf("Down sample = %d\n", dec);
    
    N = ND/dec;
    
    if(dec > 1)
    {

    fmn = 1.0;
    
    // apply 8th order zero phase bandpass filter
    bwbpf(dataF, dataF, 1, ND, 8, 1.0/dt, fmx, fmn);
    bwbpf(dataF, dataF, -1, ND, 8, 1.0/dt, fmx, fmn);
        
    }
    
    
    times = (double*)malloc(sizeof(double)* (N));
    data = (double*)malloc(sizeof(double)* (N));
    Hf = (double*)malloc(sizeof(double)* (N));
    
    // decimate
    for (i = 0; i < N; ++i)
    {
        times[i] = timeF[i*dec];
        data[i] = dataF[i*dec];
        Hf[i] = 0.0;
    }
    
    free(timeF);
    free(dataF);

    
    // reset sample rate and Nyquist
    dt = times[1]-times[0];
    fny = 1.0/(2.0*dt);
    
    if(verbose == 1)
       {
      sprintf(command, "framed_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
      out = fopen(command,"w");
      for (i = 0; i < N; ++i) fprintf(out,"%.16e %.16e\n", times[i], data[i]);
      fclose(out);
       }
    
    Ns = (int)(Tobs*fmx);
    NC = NCC+NCH;
    
     printf("max threads = %d\n", omp_get_max_threads());
    
    lf = double_vector(maxmod);
    lh = double_vector(maxmod);
    lQ = double_vector(maxmod);
    lw = double_vector(maxmod);
    dfmx = double_vector(maxmod);
    smline = double_vector(maxmod);
    fit = double_vector(maxmod);
    
    SL = double_matrix(NC,Ns);
    SM = double_matrix(NC,Ns);
    SN = double_matrix(NC,Ns);
    SMA = double_vector(Ns);
    SNA = double_vector(Ns);
    PS = double_vector(Ns);

    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    const gsl_rng_type * P;
    P = gsl_rng_default;
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));
    for(i = 0 ; i<= NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    LarrayX = double_matrix(NC,Ns);
    logLx = double_vector(NC);
    who = int_vector(NC);
    heat = double_vector(NC);
    acS = int_vector(NC);
    acL = int_vector(NC);
    cS = int_vector(NC);
    cL = int_vector(NC);
    

    for (i=0; i< NC; i++) who[i] = i;
    for (i=0; i< NCC; i++) heat[i] = 1.0;
    for (i=NCC; i< NC; i++) heat[i] = 1.1*heat[i-1];

    printf("Data volume %f seconds, %f Hz\n", Tobs, fny);
    
    specest(data, Hf, N, Ns, dt, fmx, SN[0], SM[0], PS);
    
    // FFT the cleaned data
    
    // Tukey window parameter. Flat for (1-alpha) of data
       alpha = (2.0*tuke/Tobs);
       tukey(data, alpha, N);
       gsl_fft_real_radix2_transform(data, 1, N);
       x = sqrt(Tobs)/(double)(N/2);
       for (i = 0; i < N; ++i) data[i] *= x;
    
    double *freqs;
    
    freqs = (double*)malloc(sizeof(double)*(Ns));
    for (i = 0; i < Ns; ++i)   freqs[i] = (double)(i+1)/Tobs;
    
    if(verbose == 1)
    {
        spec = fopen("data.dat","w");
        for (i = 1; i < Ns; ++i) fprintf(spec,"%e %e %e\n", freqs[i], data[i], data[N-i]);
        fclose(spec);
    }
    
    if(verbose == 1)
       {
       out = fopen("spec.dat","w");
       for (i = 0; i < Ns; ++i)
       {
           fprintf(out,"%.15e %.15e %.15e %.15e\n", freqs[i], SN[0][i], SM[0][i], PS[i]);
       }
       fclose(out);
       }
    
    flag = smoothset(Ns, SM[0], freqs, Tobs, fit, smline, &Nknot);
    if(flag == 0) return -1;
    
    ffitx = double_matrix(NC,Nknot);
    Xsline = double_matrix(NC,Nknot);
    
     for (j = 0; j < NC; ++j)
       {
         for (i = 0; i < Nknot; ++i)
         {
          Xsline[j][i] = smline[i];
          ffitx[j][i] = fit[i];
         }
      }
    
    if(verbose == 1)
    {
     out = fopen("control.dat","w");
      for (i = 0; i < Nknot; ++i)
      {
          fprintf(out,"%e %e\n", ffitx[0][i], exp(Xsline[0][i]));
      }
      fclose(out);
    }
    
    // used by the smooth linear fit
    if(itype == 1)
    {
    slopex = double_matrix(NC,Nknot);
    cx = double_matrix(NC,Nknot);
    bx = double_matrix(NC,Nknot);
    delx = double_matrix(NC,Nknot);
    }
    else
    {
     // Allocate spline
     aspline = gsl_spline_alloc(gsl_interp_akima, Nknot);
     acc = gsl_interp_accel_alloc();
    }
    
    if(itype == 1)
    {
     // set up smooth line fit
     smootset(smoothline);
     setupsline(smoothline, Nknot, Xsline[0], ffitx[0], slopex[0], bx[0], cx[0], delx[0]);
        for (j = 1; j < NC; ++j)
         {
          for (i = 0; i < Nknot; ++i)
          {
            bx[j][i] = bx[0][i];
            cx[j][i] = cx[0][i];
            delx[j][i] = delx[0][i];
            slopex[j][i] = slopex[0][i];
         }
         }
    }
    
    flag = lineset(Ns, SM[0], PS, freqs, Tobs, lf, lh, lQ, lw, dfmx, &Nlines);
     if(flag == 0) return -1;
    
    linef = double_matrix(NC,Nlines);  // central frequency
    lineh = double_matrix(NC,Nlines);  // line height
    lineQ = double_matrix(NC,Nlines); // line Q
    linew = double_matrix(NC,Nlines); // line width
    deltafmax = double_matrix(NC,Nlines); // cut-off
    
    for (j = 0; j < NC; ++j)
         {
            for (i = 0; i < Nlines; ++i)
            {
             linef[j][i] = lf[i];
             lineh[j][i] = lh[i];
             lineQ[j][i] = lQ[i];
             linew[j][i] = lw[i];
             deltafmax[j][i] = dfmx[i];
            }
         }

    
     // initialize smooth spectrum
    if(itype == 1)
    {
      sline(smoothline, 0, Ns, Nknot, freqs, SM[0], ffitx[0], Xsline[0], bx[0], cx[0], delx[0]);
    }
    else
    {
        gsl_spline_init(aspline,ffitx[0],Xsline[0],Nknot);
        SM[0][0] = exp(Xsline[0][0]);
        SM[0][Ns-1] = exp(Xsline[0][Nknot-1]);
        for (i = 1; i < Ns-1; ++i) SM[0][i] = exp(gsl_spline_eval(aspline,freqs[i],acc));
    }
    
      // initialize line spectrum
      for (i = 0; i < Ns; ++i)
      {
          f = freqs[i];
          y = 0.0;
          for (j = 0; j < Nlines; ++j) y += line(f, linef[0][j], lineh[0][j], linew[0][j], deltafmax[0][j], lineQ[0][j]);
          SL[0][i] = y;
      }
      
      
       if(verbose == 1)
       {
         out = fopen("specstart.dat","w");
         for (i = 0; i < Ns; ++i)
         {
             SN[0][i] = SM[0][i]+SL[0][i];
             fprintf(out,"%e %e %e %e\n", freqs[i], SM[0][i], SL[0][i], SN[0][i]);
         }
         fclose(out);
       }
    
         for (j = 1; j < NC; ++j)
             {
                for (i = 0; i < Ns; ++i)
                {
                 SM[j][i] = SM[0][i];
                 SN[j][i] = SN[0][i];
                 SL[j][i] = SL[0][i];
                }
             }
    
           for (i = 0; i < Ns; ++i)
            {
            LarrayX[0][i] = -(log(SN[0][i]) + PS[i]/SN[0][i]);
            }
    
           logLx[0] = 0.0;
           for (i=0; i< Ns; i++) logLx[0] += LarrayX[0][i];
    
           for (j = 1; j < NC; ++j)
               {
                 logLx[j] = logLx[0];
                 for (i = 0; i < Ns; ++i) LarrayX[j][i] = LarrayX[0][i];
               }
    
    sc = 0;
    for (i = 0; i < Ns; ++i) SNA[i] = 0.0;
    
    sma = 0;
    for (i = 0; i < Ns; ++i) SMA[i] = 0.0;
    
    
    // max likelihood
    lmax = logLx[0];
    qmax = 0;
    
    MM = 200000;   // iterations of MCMC
    if(fmx > 4000.0) MM = 400000;
    
    scount = 1;
    sacc = 0;
    for (i=0; i< NC; i++)
    {
        acS[i] = 0;
        acL[i] = 0;
        cS[i] = 0;
        cL[i] = 0;
    }

    if(verbose==1) out = fopen("schain.dat","w");
    
    itime = omp_get_wtime();

       for (mc = 1; mc < MM; ++mc)
       {

           
           if(mc < (MM/2-10000) && mc%(MM/8) == 0)  // during burnin we reset the smooth component and line component
           {
           
               // free up all the memory
               if(itype == 1)
               {
               free_double_matrix(bx,NC);
               free_double_matrix(cx,NC);
               free_double_matrix(delx,NC);
               free_double_matrix(slopex,NC);
               }
               else
               {
                  gsl_spline_free (aspline);
                  gsl_interp_accel_free (acc);
               }
               free_double_matrix(ffitx,NC);
               free_double_matrix(Xsline,NC);
               free_double_matrix(lineh,NC);
               free_double_matrix(linef,NC);
               free_double_matrix(lineQ,NC);
               free_double_matrix(linew,NC);
               free_double_matrix(deltafmax,NC);
               
               for (i = 0; i < Ns; ++i) SMA[i] /= (double)(sma);
               
               // reset the spacing for the smooth part
               flag = smoothset(Ns, SMA, freqs, Tobs, fit, smline, &Nknot);
                if(flag == 0) return -1;
               
               ffitx = double_matrix(NC,Nknot);
               Xsline = double_matrix(NC,Nknot);
               
                for (j = 0; j < NC; ++j)
                  {
                    for (i = 0; i < Nknot; ++i)
                    {
                     Xsline[j][i] = smline[i];
                     ffitx[j][i] = fit[i];
                    }
                 }
               
               // used by the smooth linear fit
               if(itype == 1)
               {
               slopex = double_matrix(NC,Nknot);
               cx = double_matrix(NC,Nknot);
               bx = double_matrix(NC,Nknot);
               delx = double_matrix(NC,Nknot);
               }
               else
               {
                // Allocate spline
                aspline = gsl_spline_alloc(gsl_interp_akima, Nknot);
                acc = gsl_interp_accel_alloc();
               }

               
               if(itype == 1)
               {
                // set up smooth line fit
                smootset(smoothline);
                setupsline(smoothline, Nknot, Xsline[0], ffitx[0], slopex[0], bx[0], cx[0], delx[0]);
                   for (j = 1; j < NC; ++j)
                    {
                     for (i = 0; i < Nknot; ++i)
                     {
                       bx[j][i] = bx[0][i];
                       cx[j][i] = cx[0][i];
                       delx[j][i] = delx[0][i];
                       slopex[j][i] = slopex[0][i];
                    }
                    }
               }
               
               // we are still using the current smooth line model here
               flag = lineset(Ns, SMA, PS, freqs, Tobs, lf, lh, lQ, lw, dfmx, &Nlines);
              if(flag == 0) return -1;
               
               linef = double_matrix(NC,Nlines);  // central frequency
               lineh = double_matrix(NC,Nlines);  // line height
               lineQ = double_matrix(NC,Nlines); // line Q
               linew = double_matrix(NC,Nlines); // line width
               deltafmax = double_matrix(NC,Nlines); // cut-off
               
               for (j = 0; j < NC; ++j)
                    {
                       for (i = 0; i < Nlines; ++i)
                       {
                        linef[j][i] = lf[i];
                        lineh[j][i] = lh[i];
                        lineQ[j][i] = lQ[i];
                        linew[j][i] = lw[i];
                        deltafmax[j][i] = dfmx[i];
                       }
                    }

               
                // initialize smooth spectrum
               if(itype == 1)
               {
                 sline(smoothline, 0, Ns, Nknot, freqs, SM[0], ffitx[0], Xsline[0], bx[0], cx[0], delx[0]);
               }
               else
               {
                   gsl_spline_init(aspline,ffitx[0],Xsline[0],Nknot);
                   SM[0][0] = exp(Xsline[0][0]);
                   SM[0][Ns-1] = exp(Xsline[0][Nknot-1]);
                   for (i = 1; i < Ns-1; ++i) SM[0][i] = exp(gsl_spline_eval(aspline,freqs[i],acc));
               }
               
                 // initialize line spectrum
                 for (i = 0; i < Ns; ++i)
                 {
                     f = freqs[i];
                     y = 0.0;
                     for (j = 0; j < Nlines; ++j) y += line(f, linef[0][j], lineh[0][j], linew[0][j], deltafmax[0][j], lineQ[0][j]);
                     SL[0][i] = y;
                 }
               
                    for (j = 1; j < NC; ++j)
                        {
                           for (i = 0; i < Ns; ++i)
                           {
                            SM[j][i] = SM[0][i];
                            SN[j][i] = SN[0][i];
                            SL[j][i] = SL[0][i];
                           }
                        }
               
                      for (i = 0; i < Ns; ++i)
                       {
                       LarrayX[0][i] = -(log(SN[0][i]) + PS[i]/SN[0][i]);
                       }
               
                      logLx[0] = 0.0;
                      for (i=0; i< Ns; i++) logLx[0] += LarrayX[0][i];
               
                      for (j = 1; j < NC; ++j)
                          {
                            logLx[j] = logLx[0];
                            for (i = 0; i < Ns; ++i) LarrayX[j][i] = LarrayX[0][i];
                          }
               
               // reset the averaging
               sma = 0;
               for (i = 0; i < Ns; ++i) SMA[i] = 0.0;
               
 
               // dynamic scaling of jumps (ok during burn-in)
               
               k = who[0];
               x = (double)acS[k]/(double)(cS[k]);
               if(x > 0.7) smsc *= 2.0;
               if(x < 0.3) smsc /= 2.0;
               
               k = who[0];
               x = (double)acL[k]/(double)(cL[k]);
               if(x > 0.7) lnsc *= 2.0;
               if(x < 0.3) lnsc /= 2.0;
               
               scount = 1;
               sacc = 0;
               for (i=0; i< NC; i++)
               {
                   acS[i] = 0;
                   acL[i] = 0;
                   cS[i] = 1;
                   cL[i] = 1;
               }
               
               
               
           }
           

           alpha = gsl_rng_uniform(r);
           
           if((NC > 1) && (alpha < 0.4))  // decide if we are doing a MCMC update of all the chains or a PT swap
           {
              // chain swap
               
               for(k=0; k < NC; k++)
               {
               alpha = (double)(NC-1)*gsl_rng_uniform(r);
               j = (int)(alpha);
               if(j > NCC-1) scount++;
               beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
               alpha = gsl_rng_uniform(r);
               if(beta > alpha)
               {
                   hold = who[j];
                   who[j] = who[j+1];
                   who[j+1] = hold;
                   if(j > NCC-1) sacc++;
               }
               }
               
           }
           else      // MCMC update
           {
           
           
           #pragma omp parallel for
           for(k=0; k < NC; k++)
           {
            int q = who[k];
            update(k, q, logLx, smoothline, heat[k], smsc, lnsc, Tobs, Ns, Nknot, Nlines, freqs, PS, LarrayX[q], SM[q], SL[q], SN[q], ffitx[q], Xsline[q],  bx[q], cx[q], slopex[q], delx[q], linef[q], lineh[q], lineQ[q], linew[q], deltafmax[q], cS, cL, acS, acL, rvec[k]);
           }
               
           }

           
          if(verbose == 1)
          {
              k = who[0];
              if(mc%1000 == 0) printf("%d %e %f %f %f\n", mc, logLx[k], (double)acS[k]/(double)(cS[k]), (double)acL[k]/(double)(cL[k]), (double)sacc/(double)(scount));
          }
           
           
           if(mc >= 1000 && mc%100 == 0)
           {
            for(k=0; k < NCC; k++)
               {
                sma++;
               for (i = 0; i < Ns; ++i) SMA[i] += SM[who[k]][i];
               }
           }
           
           
           if(mc >= MM/2 && mc%200 == 0)
            {
                
             for(k=0; k < NCC; k++)
                {
                sc++;
                for (i = 0; i < Ns; ++i) SNA[i] += SN[who[k]][i];
                }
            }
                
             /*
             if(verbose == 1)
              {
                  sprintf(command, "spec_sample_%d.dat", (mc-MM/2)/200);
                  spec = fopen(command,"w");
                  for (i = 1; i < Ns; ++i) fprintf(spec,"%e %e\n", freqs[i], SN[0][i]);
                  fclose(spec);
              }
              */
               
           
          if(verbose == 1)
          {
            if(mc%100 == 0)
            {
            fprintf(out, "%d ", mc);
            for(k=0; k < NC; k++) fprintf(out, "%e ", logLx[who[k]]);
            k = who[0];
            fprintf(out, "%f %f %f\n", (double)acS[k]/(double)(cS[k]), (double)acL[k]/(double)(cL[k]), (double)sacc/(double)(scount));
            }
          }
      
       }
       if(verbose == 1) fclose(out);
           
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("\nMCMC took %f seconds\n", exec_time);
    
    for (i = 0; i < Ns; ++i) SNA[i] /= (double)(sc);

    out = fopen("PSD.dat","w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e\n", f, SNA[i], PS[i]);
    }
    fclose(out);
    
  
    whiten(data, SNA, N);
    
    x = 0.0;
    y = 0.0;
    out = fopen("white.dat","w");
    for (i = 1; i < Ns; ++i)
    {
        f = freqs[i];
        fprintf(out,"%e %e %e\n", f, data[i], data[N-i]);
        x += data[i]+data[N-i];
        y += (data[i])*(data[i])+(data[N-i])*(data[N-i]);
    }
    fclose(out);
    x /= (double)(2*(Ns-1));
    y /= (double)(2*(Ns-1));
    
    printf("mean %f sd %f\n", x, sqrt(y-x*x));
        

      free(freqs);
    
      free(smoothline->lookup);
      if(itype == 1)
      {
      free_double_matrix(bx,NC);
      free_double_matrix(cx,NC);
      free_double_matrix(delx,NC);
      free_double_matrix(slopex,NC);
      }
      else
      {
         gsl_spline_free (aspline);
         gsl_interp_accel_free (acc);
      }
      free_double_matrix(ffitx,NC);
      free_double_matrix(Xsline,NC);
      free_double_matrix(lineh,NC);
      free_double_matrix(linef,NC);
      free_double_matrix(lineQ,NC);
      free_double_matrix(linew,NC);
      free_double_matrix(deltafmax,NC);
      free_double_matrix(SM,NC);
      free_double_matrix(SL,NC);
      free_double_matrix(LarrayX,NC);
      free_double_vector(heat);
      free_double_vector(logLx);
      free_int_vector(who);
      free_double_vector(smline);
      free_double_vector(fit);
      free(PS);
      free(times);
      free(data);
      
      
    
    return 0;

}

void update(int m, int q, double *logLx, struct Smooth *smoothline, double heat, double smsc, double lnsc, double Tobs, int Ns, int Nknot, int Nlines, double *freqs, double *PS, double *LarrayX, double *SM, double *SL, double *SN, double *ffitx, double *Xsline,  double *bx, double *cx, double *slopex, double *delx,  double  *linef, double *lineh, double *lineQ, double *linew, double *deltafmax, int *cS, int *cL, int *acS, int *acL, gsl_rng *r)
{
    double alpha, H, logLy;
    int i, j, k;
    int typ;
    int imin, imax;
    double *Ysline, *ffity;
    double *slopey, *by, *cy, *dely;
    double *DS, *SNY, *LarrayY;
    double ylinew, ylinef, ylineh, ylineQ, ydeltafmax;
    double sqht, df, spread;
    
    gsl_spline   *aspline;
    gsl_interp_accel *acc;
    
    sqht = sqrt(heat);
    df = 1.0/Tobs;
    
    alpha = gsl_rng_uniform(r);
    
    SNY = double_vector(Ns);
    LarrayY = double_vector(Ns);
    DS = double_vector(Ns);
    
    for (i = 0; i < Ns; ++i)
    {
        SNY[i] = SN[i];
        LarrayY[i] = LarrayX[i];
    }
    

  if(alpha < 0.6) // update sline
  {
    
   typ = 0;
   cS[m]++;
      
      Ysline = double_vector(Nknot);
      ffity = double_vector(Nknot);
      
      if(itype == 1)
      {
          slopey = double_vector(Nknot);
          by = double_vector(Nknot);
          cy = double_vector(Nknot);
          dely = double_vector(Nknot);
      }
      else
      {
        // Allocate spline
        aspline = gsl_spline_alloc(gsl_interp_akima, Nknot);
        acc = gsl_interp_accel_alloc();
      }
      
   
     // important to make sure all the sline points and slopes are properly initialized
     // the proposed updates reset parts of slopey array, and the impact on the aa, bb etc
     // coefficients extends over several points
     for (i = 0; i < Nknot; ++i)
     {
         Ysline[i] =  Xsline[i];
         ffity[i] = ffitx[i]; // the knot locations don't change in this code
     }
        
     if(itype == 1)
     {
       for (i = 0; i < Nknot; ++i)
         {
          slopey[i] = slopex[i];
          by[i] = bx[i];
          cy[i] = cx[i];
          dely[i] = delx[i];
         }
     }
        
      // pick a knot to update
      k = (int)((double)(Nknot)*gsl_rng_uniform(r));
        
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.7)
        {
           Ysline[k] =  Xsline[k] + smsc*sqht*gsl_ran_gaussian(r,0.1);
        }
        else if (alpha > 0.4)
        {
            Ysline[k] =  Xsline[k] + smsc*sqht*gsl_ran_gaussian(r,0.5);
        }
        else
        {
            Ysline[k] =  Xsline[k] + smsc*sqht*gsl_ran_gaussian(r,2.0);
        }
        

        if(itype == 1)
        {
        // if control point k is updated, need to do a delta update on the smooth line model
        // region between ffit[k-(Nsten-2)] and ffit[k] is impacted. Care needs to be taken at boundary points.
        delta_setupsline(smoothline, k, Nknot, Ysline, ffity, slopey, by, cy, dely);
        getrangesmooth(smoothline, k, Nknot, Ns, ffity, Tobs, &imin, &imax);
        sline(smoothline, imin, imax, Nknot, freqs, DS, ffity, Ysline, by, cy, dely);
        }
        else
        {
          gsl_spline_init(aspline,ffitx,Ysline,Nknot);
          getrangeakima(k, Nknot, ffity, Tobs, &imin, &imax);
          for (i = imin; i < imax; ++i) DS[i] = exp(gsl_spline_eval(aspline,freqs[i],acc));
        }
        
       
        for (i = imin; i < imax; ++i)
         {
           SNY[i] = DS[i] + SL[i];
           LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
         }
        logLy = logLx[q];
        for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
      
      
         
        
    }
    else // line update
    {
        
       typ = 1;
       cL[m]++;
        
        // pick a line to update
                   k = (int)((double)(Nlines)*gsl_rng_uniform(r));
                    
                    alpha = gsl_rng_uniform(r);
                    if(alpha > 0.5)
                    {
                     ylinef = linef[k] + lnsc*sqht*gsl_ran_gaussian(r,df);
                     ylineh = lineh[k]*(1.0+lnsc*sqht*gsl_ran_gaussian(r,0.05));
                     ylineQ = lineQ[k] + lnsc*sqht*gsl_ran_gaussian(r,1.0);
                    }
                    else if (alpha > 0.3)
                    {
                     ylinef = linef[k] + lnsc*sqht*gsl_ran_gaussian(r,2.0*df);
                     ylineh = lineh[k]*(1.0+lnsc*sqht*gsl_ran_gaussian(r,0.1));
                     ylineQ = lineQ[k] + lnsc*sqht*gsl_ran_gaussian(r,2.0);
                    }
                    else
                    {
                     ylinef = linef[k] + lnsc*sqht*gsl_ran_gaussian(r,0.2*df);
                     ylineh = lineh[k]*(1.0+lnsc*sqht*gsl_ran_gaussian(r,0.01));
                     ylineQ = lineQ[k] + lnsc*sqht*gsl_ran_gaussian(r,0.2);
                    }
        
       if(ylineQ < 10.0) ylineQ = 10.0;
       
        spread = (1.0e-2*ylineQ);
        if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/20
        ydeltafmax = ylinef/spread;
        ylinew = 8.0*ydeltafmax;
        
       // need to cover old and new line
        imin = (int)((linef[k]-linew[k])*Tobs);
        imax = (int)((linef[k]+linew[k])*Tobs);
        i = (int)((ylinef-ylinew)*Tobs);
        if(i < imin) imin = i;
        i = (int)((ylinef+ylinew)*Tobs);
        if(i > imax) imax = i;
        if(imin < 0) imin = 0;
        if(imax > Ns-1) imax = Ns-1;
        
       // proposed line contribution
        for (i = imin; i <= imax; ++i) DS[i] = line(freqs[i], ylinef, ylineh, ylinew, ydeltafmax, ylineQ);
          
       // have to recompute and remove current line since lines can overlap
        for (i = imin; i <= imax; ++i) DS[i] -= line(freqs[i], linef[k], lineh[k], linew[k], deltafmax[k], lineQ[k]);
        
       
        for (i = imin; i < imax; ++i)
         {
           SNY[i] = SM[i] + SL[i] + DS[i];
           LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
         }
        logLy = logLx[q];
        for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
        
    }


    H = (logLy-logLx[q])/heat;
    
    alpha = log(gsl_rng_uniform(r));

    if(H > alpha)
    {
        logLx[q] = logLy;
        
        if(typ == 0)
        {
            
        acS[m]++;
        Xsline[k] = Ysline[k];
        ffitx[k] = ffity[k];
            
        for (i = imin; i < imax; ++i)
           {
               SM[i] = DS[i]; // replace
               SN[i] = SNY[i];
               LarrayX[i] = LarrayY[i];
           }
         
            if(itype == 1)
            {
            for(i=0; i< Nknot; i++)
             {
               slopex[i] = slopey[i];
               cx[i] = cy[i];
               bx[i] = by[i];
               delx[i] = dely[i];
             }
            }
            
    
         free_double_vector(Ysline);
         free_double_vector(ffity);
         
         if(itype == 1)
         {
             free_double_vector(slopey);
             free_double_vector(by);
             free_double_vector(cy);
             free_double_vector(dely);
         }
         else
         {
           // free spline
              gsl_spline_free (aspline);
              gsl_interp_accel_free (acc);
         }
            
        }
        
        if(typ == 1)
        {
        acL[m]++;
          
           linef[k] = ylinef;
           lineh[k] = ylineh;
           lineQ[k] = ylineQ;
           linew[k] = ylinew;
           deltafmax[k] = ydeltafmax;
           
          for (i = imin; i <= imax; ++i)
           {
               SL[i] += DS[i];   // delta this segment
               SN[i] = SNY[i];
               LarrayX[i] = LarrayY[i];
           }
            
        }
        
    }
    
    free_double_vector(DS);
    free_double_vector(SNY);
    free_double_vector(LarrayY);
    

}

int smoothset(int Ns, double *SM, double *freqs, double Tobs, double *fit, double *line, int *Nk)
{
    int i, sm;
    int j, k, ii, kk, flag;
    double x, max;
    int Nknot;
    double *S1, *S2;
    
    S1 = (double*)malloc(sizeof(double)*(Ns));
    S2 = (double*)malloc(sizeof(double)*(Ns));
    
     // moving average
     sm = (int)(smooth*Tobs);
     x = 0.0;
     for (i = 0; i < sm; ++i) x += SM[i];
     for (i = sm; i < Ns; ++i)
     {
         S1[i-sm/2] = x/(double)(sm);
         x += SM[i] - SM[i-sm];
     }
     
     // moving average with wider window
     sm *= 2;
     x = 0.0;
     for (i = 0; i < sm; ++i) x += SM[i];
     for (i = sm; i < Ns; ++i)
     {
         S2[i-sm/2] = x/(double)(sm);
         x += SM[i] - SM[i-sm];
     }
     
     // fill initial bins
     for (i = 0; i < sm/2; ++i)
     {
         S1[i] = SM[i];
         S2[i] = 2.0*S1[i];
     }
     
    // count the number of spline knots
     Nknot = 1;
     k = (int)(dfmin*Tobs);
     kk = (int)(dfmax*Tobs);
     j = 0;
     flag = 0;
     max = 0.0;
     for (i = 1; i < Ns; ++i)
     {
         x = fabs(S2[i]/S1[i]-1.0);
         if(x > max) max = x;
         j++;
          if(i%k == 0)
           {
               if(max > tol || j == kk)
               {
                  // printf("%f %f %f\n", (double)(i-j/2)/Tobs, (double)(j)/Tobs, max);
                   max = 0.0;
                   j = 0;
                   Nknot++;
               }
           }
         
         
     }
     
      Nknot++;
     
      printf("There are %d spline knots\n", Nknot);
    
      if(Nknot > maxmod)
      {
          printf("To many points in the smooth fit\n");
          return 0;
      }
      
      fit[0] = freqs[0];
      line[0] = log(S1[0]);
      
      ii = 1;
      j = 0;
      flag = 0;
      max = 0.0;
      for (i = 1; i < Ns; ++i)
      {
          x = fabs(S2[i]/S1[i]-1.0);
          if(x > max) max = x;
          j++;
           if(i%k == 0)
            {
                if(max > tol || j == kk)
                {
                    max = 0.0;
                    fit[ii] = freqs[(i-j/2)];
                    line[ii] = log(S1[i-j/2]);
                    j = 0;
                    ii++;
                }
            }
          
          
      }
      fit[Nknot-1] = freqs[Ns-1];
      line[Nknot-1] = log(SM[Ns-1]);

    *Nk = Nknot;
    
    free_double_vector(S1);
    free_double_vector(S2);
    
    return 1;
    
}

int lineset(int Ns, double *SM, double *PS, double *freqs, double Tobs, double *lf, double *lh, double *lQ, double *lw, double *dfmx, int *Nlns)
{
    int Nlines, flag, i, j, k, ii;
    double x, max, xold, spread;
    
    // count the number of lines
   j = 0;
   flag = 0;
   for (i = 0; i < Ns; ++i)
   {
       x = PS[i]/SM[i];
       // start of a line
       if(x > linemul && flag == 0)
       {
           k = 1;
           flag = 1;
           max = x;
           ii = i;
       }
       // in a line
       if(x > linemul  && flag ==1)
       {
           k++;
           if(x > max)
           {
               max = x;
               ii = i;
           }
       }
       // have reached the end of a line
       if(flag == 1)
       {
           if(x < linemul)
           {
               flag = 0;
               j++;
           }
       }
   }
  
   
   Nlines = j;
    
    if(Nlines > maxmod)
    {
        printf("To many lines in the line model\n");
        return 0;
    }
  
  j = -1;
  xold = 1.0;
  flag = 0;
  for (i = 0; i < Ns; ++i)
  {
      x = PS[i]/SM[i];
      // start of a line
      if(x > linemul && flag == 0)
      {
          k = 1;
          flag = 1;
          max = x;
          ii = i;
      }
      // in a line
      if((x > linemul) && flag ==1)
      {
          k++;
          if(x > max)
          {
              max = x;
              ii = i;
          }
      }
      // have reached the end of a line
      if(flag == 1)
      {
          if(x < linemul)
          {
              flag = 0;
              j++;
              lf[j] = freqs[ii];
              lh[j] = (max-1.0)*SM[ii];
              lQ[j] = sqrt(max)*lf[j]*Tobs/(double)(k);
              
                spread = (1.0e-2*lQ[j]);
                if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/50
                dfmx[j] = lf[j]/spread;
                lw[j] = 8.0*dfmx[j];
             
          }
      }

  }
  
  printf("There are %d lines\n", Nlines);
    
    *Nlns = Nlines;
    
    return 1;

}
