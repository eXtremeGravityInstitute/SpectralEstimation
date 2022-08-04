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

int main(int argc, char *argv[])
{

  int i, j, k, ii, kk, sc, ND, N, Ns, Nm, m, mc, dec, MM, MK, MS, skip;
  int scount, sacc, hold;
  int NC, NCC, NCH;
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

    
    // omp can count threads, but doesn't seem to know about physical cores
    // hyperthreading does appear to help MPI, so physical cores are what really count
    // Unfortunately, OSX have different and incompatible system calles to tell you
    // the real number of cores (sysctl -n hw.physicalcpu in OSX and lscpu in Linux)
    // Instead, we divide down by 2 to be safe
    NC = omp_get_max_threads()/2;
    if(NC > coremax) NC = coremax; // don't wnat to be too greedy
    NCC = NC/2;  // use half the cores for cold chains
    NCH = NC-NCC;
    
    printf("Using %d cold and %d hot chains\n", NCC, NCH);
    
    
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
             if(addmul == 0)
             {
                 SN[0][i] = SM[0][i]+SL[0][i];
             }
             else
             {
                 SN[0][i] = SM[0][i]*(1.0+SL[0][i]);
             }
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
    
    MM = 200000/NCC;   // iterations of MCMC
    if(fmx > 4000.0) MM = 400000/NCC;
    
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
                     if(addmul == 0)
                       {
                        SN[0][i] = SM[0][i]+SL[0][i];
                       }
                      else
                       {
                        SN[0][i] = SM[0][i]*(1.0+SL[0][i]);
                       }
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
               if(x < 0.4) smsc /= 2.0;
               
               k = who[0];
               x = (double)acL[k]/(double)(cL[k]);
               if(x > 0.7) lnsc *= 2.0;
               if(x < 0.4) lnsc /= 2.0;
               
               printf("%f %f\n", smsc, lnsc);
               
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

