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


//OSX
// clang -O3 -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c SpecUtils.c -lgsl  -lm

// Linux
// gcc -O3 -std=gnu99 -fopenmp -w -o SpecFit SpecFit.c SpecUtils.c -lgsl -lm

int main(int argc, char *argv[])
{

  int i, j, k, ii, kk, sc, ND, N, Ns, Nm, m, mc, dec, MK, MS, skip;
  int imin, imax, acS, acL, acH, cS, cL, cH;
  int typ, sm;
  double *timeF, *dataF;
  double *times, *data, *Hf;
  double x, y, z, dt, df, Tobs, ttrig, fny, fmn, fmx;
  double f, finc;
  double *SN, *SNY, *SNA, *SM, *SL, *PS, *S1, *S2, *sp;
  double *LarrayX, *LarrayY;
  int Nsp, Nl;
  char command[1024];
  double spread, *deltafmax, *linew;
  double *linef, *lineh, *lineQ;
  double f2, f4, ff4, ff2, df2, df4;
  double ylinew, ylinef, ylineh, ylineQ, ydeltafmax;
  int flag, Nlines;
  double xold, xnext, Abar;
  double max;
  int Nknot, Nsten;
  double a1, a2, sdd;
    
  // these are used by smooth linear fit. Have to be declared even if not allocated
  double *bx, *cx, *delx;
  double *by, *cy, *dely;
  double *slopex, *slopey;
  struct Smooth *smoothline  = malloc(sizeof(struct Smooth));
    
  // these are used by the spline. Have to be declared even if not allocated
  gsl_spline   *aspline;
  gsl_interp_accel *acc;
    
  double alpha;
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
    
    printf("max threads = %d\n", omp_get_max_threads());
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NTHREAD);
    //##############################################
    
    
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

    fmn = 8.0;
    
    // apply 8th order zero phase bandpass filter
    bwbpf(dataF, dataF, 1, ND, 8, 1.0/dt, fmx, fmn);
    bwbpf(dataF, dataF, -1, ND, 8, 1.0/dt, fmx, fmn);
    
    N = ND/dec;
    
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
    
    S1 = (double*)malloc(sizeof(double)*(Ns));
    S2 = (double*)malloc(sizeof(double)*(Ns));
    SL = (double*)malloc(sizeof(double)*(Ns));
    SM = (double*)malloc(sizeof(double)*(Ns));
    SN = (double*)malloc(sizeof(double)*(Ns));
    SNA = (double*)malloc(sizeof(double)*(Ns));
    SNY = (double*)malloc(sizeof(double)*(Ns));
    PS = (double*)malloc(sizeof(double)*(Ns));
    LarrayX = (double*)malloc(sizeof(double)*(Ns));
    LarrayY = (double*)malloc(sizeof(double)*(Ns));
    
    printf("Data volume %f seconds, %f Hz\n", Tobs, fny);
    
    specest(data, Hf, N, Ns, dt, fmx, SN, SM, PS);
    
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
           fprintf(out,"%.15e %.15e %.15e %.15e\n", freqs[i], SN[i], SM[i], PS[i]);
       }
       fclose(out);
       }
    
    
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

    /*
    out = fopen("smooth.dat","w");
    for (i = sm/2; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e\n", f, S1[i], S2[i]);
    }
    fclose(out);
    */
    
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
    
     double *ffitx, *ffity, *Xsline, *Ysline;
     
     ffitx = (double*)malloc(sizeof(double)*(Nknot));
     ffity = (double*)malloc(sizeof(double)*(Nknot));
     Xsline = (double*)malloc(sizeof(double)*(Nknot));
     Ysline = (double*)malloc(sizeof(double)*(Nknot));
    
    // used by the smooth linear fit
    if(itype == 1)
    {

    slopex = (double*)malloc(sizeof(double)*(Nknot));
    cx = (double*)malloc(sizeof(double)*(Nknot));
    bx = (double*)malloc(sizeof(double)*(Nknot));
    delx = (double*)malloc(sizeof(double)*(Nknot));
    
    slopey = (double*)malloc(sizeof(double)*(Nknot));
    cy = (double*)malloc(sizeof(double)*(Nknot));
    by = (double*)malloc(sizeof(double)*(Nknot));
    dely = (double*)malloc(sizeof(double)*(Nknot));
    }
    
    ffitx[0] = freqs[0];
    Xsline[0] = log(S1[0]);
    
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
                  ffitx[ii] = freqs[(i-j/2)];
                  Xsline[ii] = log(S1[i-j/2]);
                  j = 0;
                  ii++;
              }
          }
        
        
    }
    ffitx[Nknot-1] = freqs[Ns-1];
    Xsline[Nknot-1] = log(SM[Ns-1]);
    
    if(verbose == 1)
    {
     out = fopen("control.dat","w");
      for (i = 0; i < Nknot; ++i)
      {
          fprintf(out,"%e %e\n", ffitx[i], exp(Xsline[i]));
      }
      fclose(out);
    }

    
    if(itype == 1)
    {
     // set up smooth line fit
     smootset(smoothline);
     setupsline(smoothline, Nknot, Xsline, ffitx, slopex, bx, cx, delx);
      for (i = 0; i < Nknot; ++i)
       {
            by[i] = bx[i];
            cy[i] = cx[i];
            dely[i] = delx[i];
            slopey[i] = slopex[i];
       }
    }
    else
    {
     // Allocate spline
     aspline = gsl_spline_alloc(gsl_interp_akima, Nknot);
     acc = gsl_interp_accel_alloc();
     /* compute spline */
     gsl_spline_init(aspline,ffitx,Xsline,Nknot);
    }

   

     for (i = 0; i < Nknot; ++i)
        {
            Ysline[i] = Xsline[i];
            ffity[i] = ffitx[i];  // These never change (we don't update the locations or number of control points in this code)
        }
    
      

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
      


         linef = (double*)malloc(sizeof(double)*(Nlines));  // central frequency
         lineh = (double*)malloc(sizeof(double)*(Nlines));  // line height
         lineQ = (double*)malloc(sizeof(double)*(Nlines)); // line Q
         linew = (double*)malloc(sizeof(double)*(Nlines));  // line width
         deltafmax = (double*)malloc(sizeof(double)*(Nlines));  // cut-off
      
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
                  linef[j] = freqs[ii];
                  lineh[j] = (max-1.0)*SM[ii];
                  lineQ[j] = sqrt(max)*linef[j]*Tobs/(double)(k);
                  
                    spread = (1.0e-2*lineQ[j]);
                    if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/50
                    deltafmax[j] = linef[j]/spread;
                    linew[j] = 8.0*deltafmax[j];
                  
                  //printf("%d %e %e %e %e %e\n", j, linef[j], lineh[j], linew, freqs[ii], lineQ[j]);
                 
              }
          }

      }
      
      printf("There are %d lines\n", Nlines);
    
     // initialize smooth spectrum
    if(itype == 1)
    {
      sline(smoothline, 0, Ns, Nknot, freqs, SM, ffitx, Xsline, bx, cx, delx);
    }
    else
    {
        gsl_spline_init(aspline,ffitx,Xsline,Nknot);
        SM[0] = exp(Xsline[0]);
        SM[Ns-1] = exp(Xsline[Nknot-1]);
        for (i = 1; i < Ns-1; ++i) SM[i] = exp(gsl_spline_eval(aspline,freqs[i],acc));
    }
    
      // initialize line spectrum
      //out = fopen("lines.dat","w");
      for (i = 0; i < Ns; ++i)
      {
          f = freqs[i];
          y = 0.0;
          for (j = 0; j < Nlines; ++j) y += line(f, linef[j], lineh[j], linew[j], deltafmax[j], lineQ[j]);
          SL[i] = y;
          //fprintf(out,"%e %e %e\n", f, SM[i] + SL[i], SL[i]);
      }
     //fclose(out);
      
      
       if(verbose == 1)
       {
         out = fopen("specstart.dat","w");
         for (i = 0; i < Ns; ++i)
         {
             fprintf(out,"%e %e %e %e\n", freqs[i], SM[i], SL[i], SM[i]+SL[i]);
         }
         fclose(out);
       }
    
      // This holds the updated values
      
      double *DSM, *DSL;
      
       DSM = (double*)malloc(sizeof(double)*(Ns));
       DSL = (double*)malloc(sizeof(double)*(Ns));
     
      double logLx, logLy, logpx, logpy, pmulx, pmuly;
      double H;
    
     // initialize delta arrays
        for (i = 0; i < Ns; ++i)
         {
          DSM[i] = SM[i];
          DSL[i] = 0.0;
        }
    
        if(ompflag == 1)
            {
                // initialize likelihood
                #pragma omp parallel for
                for (i = 0; i < Ns; ++i)
                {
                  SN[i] = SM[i] + SL[i];
                  SNY[i] = SN[i];
                  LarrayX[i] = -(log(SN[i]) + PS[i]/SN[i]);
                }
                
                logLx = 0.0;
                for (i=0; i< Ns; i++) logLx += LarrayX[i];

                /*
                logLx = 0.0;
                #pragma omp parallel for reduction (+:logLx)
                 {
                     for (i=0; i< Ns; i++) logLx += LarrayX[i];
                 } */
                
               }
                else
              {
                  
                // non omp version
                for (i = 0; i < Ns; ++i)
                  {
                     SN[i] = SM[i] + SL[i];
                     SNY[i] = SN[i];
                     LarrayX[i] = -(log(SN[i]) + PS[i]/SN[i]);
                   }
                 logLx = 0.0;
                 for (i=0; i< Ns; i++) logLx += LarrayX[i];
                }
      
      acS = 0;
      acL = 0;
      acH = 0;
      cS = 1;
      cL = 1;
      cH = 1;
    
    double a, b;
    
    a = 0.5;
    b = 1.0;
    
    sc = 0.0;
    for (i = 0; i < Ns; ++i) SNA[i] = 0.0;

    if(verbose==1) out = fopen("schain.dat","w");
    
    itime = omp_get_wtime();
      
       for (mc = 0; mc < MM; ++mc)
       {
           
           // prune any weak lines
           if(mc > MM/4 && mc < MM/2 && mc%2000 == 0)
           {
               j = 0;
               for (k = 0; k < Nlines; ++k)
               {
                   i = (int)(linef[k]*Tobs);  // bin where line peaks
                   x = lineh[k]/SM[i];     // height of line relative to smooth
                   if(x > 3.0)  // keep this line
                   {
                       linef[j] = linef[k];
                       lineh[j] = lineh[k];
                       lineQ[j] = lineQ[k];
                       linew[j] = linew[k];
                       deltafmax[j] = deltafmax[k];
                       j++;
                   }
               }
               
               Nlines = j;
               
               // reset line spectrum
               for (i = 0; i < Ns; ++i)
               {
                  f = freqs[i];
                  y = 0.0;
                  for (j = 0; j < Nlines; ++j) y += line(f, linef[j], lineh[j], linew[j], deltafmax[j], lineQ[j]);
                  SL[i] = y;
               }
               
                 if(ompflag == 1)
                 {
                 // reset likelihood
                 #pragma omp parallel for
                 for (i = 0; i < Ns; ++i)
                  {
                      SN[i] = SM[i] + SL[i];
                      SNY[i] = SN[i];
                      LarrayX[i] = -(log(SN[i]) + PS[i]/SN[i]);
                  }
                     
                     logLx = 0.0;
                     for (i=0; i< Ns; i++) logLx += LarrayX[i];
                     
                     /*
                    logLx = 0.0;
                    #pragma omp parallel for reduction (+:logLx)
                    {
                     for (i=0; i< Ns; i++) logLx += LarrayX[i];
                    }*/

                 }
                 else
                 {
               
               // non omp version
                 for (i = 0; i < Ns; ++i)
                  {
                      SN[i] = SM[i] + SL[i];
                      SNY[i] = SN[i];
                      LarrayX[i] = -(log(SN[i]) + PS[i]/SN[i]);
                  }
                 logLx = 0.0;
                 for (i=0; i< Ns; i++) logLx += LarrayX[i];
                 }

               
            }
             
           
           alpha = gsl_rng_uniform(r);
           
           if(alpha < a) // update sline
           {
               
            typ = 0;
            cS++;
              
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
               if(alpha > 0.5)
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.05);
               }
               else if (alpha > 0.3)
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.02);
               }
               else
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.1);
               }
               

               if(itype == 1)
               {
               // if control point k is updated, need to do a delta update on the smooth line model
               // region between ffit[k-(Nsten-2)] and ffit[k] is impacted. Care needs to be taken at boundary points.
               delta_setupsline(smoothline, k, Nknot, Ysline, ffity, slopey, by, cy, dely);
               getrangesmooth(smoothline, k, Nknot, Ns, ffity, Tobs, &imin, &imax);
               sline(smoothline, imin, imax, Nknot, freqs, DSM, ffity, Ysline, by, cy, dely);
               }
               else
               {
                 gsl_spline_init(aspline,ffitx,Ysline,Nknot);
                 getrangeakima(k, Nknot, ffity, Tobs, &imin, &imax);
                 for (i = imin; i < imax; ++i) DSM[i] = exp(gsl_spline_eval(aspline,freqs[i],acc));
               }
               
               if(ompflag == 1)
                   {
                 // compute change in likelihood
                #pragma omp parallel for
                for (i = imin; i < imax; ++i)
                 {
                  SNY[i] = DSM[i] + SL[i];
                  LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
                 }
                       
                logLy = logLx;
                for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
                       
                  /*
                 logLy = logLx;
                #pragma omp parallel for reduction (+:logLy)
                 {
                  for (i=imin; i< imax; i++) {logLy += (LarrayY[i]-LarrayX[i]);}
                 }*/
                       
                   }
                   else
                   {
               
               // non omp version
               for (i = imin; i < imax; ++i)
                {
                  SNY[i] = DSM[i] + SL[i];
                  LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
                }
               logLy = logLx;
               for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
                   }
                
               
           }
           else // line update
           {
               
              typ = 1;
              cL++;
               
               // pick a line to update
                          k = (int)((double)(Nlines)*gsl_rng_uniform(r));
                           
                           alpha = gsl_rng_uniform(r);
                           if(alpha > 0.5)
                           {
                            ylinef = linef[k] + gsl_ran_gaussian(r,df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.05));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,1.0);
                           }
                           else if (alpha > 0.3)
                           {
                            ylinef = linef[k] + gsl_ran_gaussian(r,4.0*df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.2));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,4.0);
                           }
                           else
                           {
                            ylinef = linef[k] + gsl_ran_gaussian(r,0.2*df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.01));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,0.2);
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
               for (i = imin; i <= imax; ++i) DSL[i] = line(freqs[i], ylinef, ylineh, ylinew, ydeltafmax, ylineQ);
                 
              // have to recompute and remove current line since lines can overlap
               for (i = imin; i <= imax; ++i) DSL[i] -= line(freqs[i], linef[k], lineh[k], linew[k], deltafmax[k], lineQ[k]);
               
               // compute change in likelihood
               if(ompflag == 1)
                   {
               #pragma omp parallel for
               for (i = imin; i < imax; ++i)
                 {
                   SNY[i] = SM[i] + SL[i] + DSL[i];
                   LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
                 }
            
                /*
                logLy = logLx;
                 #pragma omp parallel for reduction (+:logLy)
                 {
                    for (i=imin; i< imax; i++) {logLy += (LarrayY[i]-LarrayX[i]);}
                 } */
                       
                       logLy = logLx;
                       for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
                
                   }
                  else
                  {
               
               // non omp version
               for (i = imin; i < imax; ++i)
                {
                  SNY[i] = SM[i] + SL[i] + DSL[i];
                  LarrayY[i] =  -(log(SNY[i]) + PS[i]/SNY[i]);
                }
               logLy = logLx;
               for (i=imin; i< imax; i++) logLy += (LarrayY[i]-LarrayX[i]);
                  }
               
               
           }
           
           
           H = (logLy-logLx);
           alpha = log(gsl_rng_uniform(r));
           
           if(H > alpha)
           {
               logLx = logLy;
               
               if(typ == 0)
               {
                   
               acS++;
               Xsline[k] = Ysline[k];
               ffitx[k] = ffity[k];
                   
               for (i = imin; i < imax; ++i)
                  {
                      SM[i] = DSM[i]; // replace
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
                   
               }
               
               if(typ == 1)
               {
               acL++;
                 
                  linef[k] = ylinef;
                  lineh[k] = ylineh;
                  lineQ[k] = ylineQ;
                  linew[k] = ylinew;
                  deltafmax[k] = ydeltafmax;
                  
                 for (i = imin; i <= imax; ++i)
                  {
                      SL[i] += DSL[i];   // delta this segment
                      SN[i] = SNY[i];
                      LarrayX[i] = LarrayY[i];
                  }
                   
               }
                 
           }
           
          if(verbose == 1)  if(mc%1000 == 0) printf("%d %e %f %f\n", mc, logLx, (double)acS/(double)(cS), (double)acL/(double)(cL));
           
           if(mc >= MM/2 && mc%200 == 0)
            {
              sc++;
              for (i = 0; i < Ns; ++i) SNA[i] += SN[i];
                
             if(verbose == 1)
              {
                  sprintf(command, "spec_sample_%d.dat", (mc-MM/2)/200);
                  spec = fopen(command,"w");
                  for (i = 1; i < Ns; ++i) fprintf(spec,"%e %e\n", freqs[i], SN[i]);
                  fclose(spec);
              }
                
            }
               
           
          if(verbose == 1)  if(mc%100 == 0) fprintf(out, "%d %e %f %f\n", mc, logLx, (double)acS/(double)(cS), (double)acL/(double)(cL));
      
       }
       if(verbose == 1) fclose(out);
           
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("\nMCMC took %f seconds\n", exec_time);
        
      printf("There are %d lines remaining\n", Nlines);
    
    for (i = 0; i < Ns; ++i) SNA[i] /= (double)(sc);

    out = fopen("PSD.dat","w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e\n", f, SNA[i], PS[i]);
    }
    fclose(out);
    
  
    whiten(data, SN, N);
    
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
      free(bx);
      free(cx);
      free(delx);
      free(slopex);
      free(by);
      free(cy);
      free(dely);
      free(slopey);
      free(DSM);
      free(DSL);
      free(lineh);
      free(linef);
      free(lineQ);
      free(linew);
      free(deltafmax);
      free(SM);
      free(SN);
      free(SNY);
      free(PS);
      free(LarrayX);
      free(LarrayY);
      free(times);
      free(data);
      
      
    
    return 0;

}


