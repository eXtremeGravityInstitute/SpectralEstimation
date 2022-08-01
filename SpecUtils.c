/*******************************************************************************************

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

**********************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "SpecUtils.h"
#include "Constants.h"

void smootset(struct Smooth *smoothline)
{
    int i;
    double x, xx;
    // asmooth min for Nsten = 3 is 15
    // asmooth min for Nsten = 5 is 8
    // Note that Nknot must be >= Nsten
    // recommended settings ar Nsten =5 and asmooth between 8 and 15
    smoothline->Nsten = 5;
    smoothline->Nlook = 2000;
    smoothline->asmooth = 10.0;
    smoothline->lookup = (double*)malloc(sizeof(double)*(smoothline->Nlook));
    xx = 6.0/(double)(smoothline->Nlook);
    smoothline->dx = xx;
    for(i=0; i< smoothline->Nlook; i++)
    {
        x = (double)(i)*xx;
        smoothline->lookup[i] = log(2.0*cosh(x));
    }
}

void getrangesmooth(struct Smooth *smoothline, int k, int Nknot, int Ns, double *ffit, double Tobs, int *imin, int *imax)
{
              if(k >= (smoothline->Nsten-2))
              {
              *imin = (int)(ffit[k-(smoothline->Nsten-2)]*Tobs);
              }
              else
              {
              *imin = 0;
              }
              
              if(k < Nknot-1)
               {
               *imax = (int)(ffit[k+1]*Tobs)+1;
               }
              else
               {
               *imax = Ns;
               }
               
               // start is special
               if(k < smoothline->Nsten) *imax = (int)(ffit[smoothline->Nsten]*Tobs)+1;
    
              if(*imin < 0) *imin = 0;
              if(*imax > Ns) *imax = Ns;
}

void getrangeakima(int k, int Nknot, double *ffit, double Tobs, int *imin, int *imax)
{
              if(k > 1)
              {
                  *imin = (int)(ffit[k-2]*Tobs);
              }
              else
              {
                  *imin = (int)(ffit[0]*Tobs);
              }
              
              if(k < Nknot-2)
              {
                  *imax = (int)(ffit[k+2]*Tobs);
              }
              else
              {
                  *imax = (int)(ffit[Nknot-1]*Tobs);
              }
}

double line(double f, double linef, double lineh, double linew, double deltafmax, double lineQ)
{
    double x, y, z;
    double ff2, f2, df2, df4;
    
             f2 = linef*linef;
    
             y = 0.0;
             x = fabs(f - linef);
             if(x < linew)
             {
             z = 1.0;
             if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
             ff2 = f*f;
             df2 = lineQ*(1.0-ff2/f2);
             df4 = df2*df2;
             y = z*lineh/(ff2/f2+df4);
             }
    
    return y;
                 
}

double ltc(struct Smooth *smoothline, double x)
{
    double y;
    int i;
    
    i = (int)(fabs(x)/smoothline->dx);
    
    if(i > smoothline->Nlook-1)  // ouside of lookup
    {
        y = log(2.0*cosh(x));
    }
    else
    {
        y = smoothline->lookup[i];
    }
    
    return y;
    
}


void setupsline(struct Smooth *smoothline, int Nknot, double *sline, double *ffit, double *slopes, double *bb, double *cc, double *dd)
{
    int i, k;
    
      for(i=0; i< Nknot-1; i++)
        {
         slopes[i] = (sline[i+1]-sline[i])/(ffit[i+1]-ffit[i]);
         dd[i] = 0.5/(ffit[i+1]-ffit[i]);
        }
    
      for(i=0; i< Nknot-2; i++)
       {
           cc[i+1] = 0.5*(slopes[i+1]-slopes[i]);
           k = i+(smoothline->Nsten-2);
           if(k > Nknot-2) k = Nknot-2;
           bb[i] = 0.5*(slopes[k]+slopes[i]);
       }
        
}

// when point iu is updated, the fit between control points at iu-(Nsten-1)/2;  and iu+(Nsten-1)/2 is changed
void delta_setupsline(struct Smooth *smoothline, int iu, int Nknot, double *sline, double *ffit, double *slopes,  double *bb, double *cc, double *dd)
{
    int ii, i, k;

    
      for(i=-1; i< 1; i++)
       {
        ii = iu+i;
        if(ii > -1 && ii < Nknot-1)
        {
        slopes[ii] = (sline[ii+1]-sline[ii])/(ffit[ii+1]-ffit[ii]);
        dd[ii] = 0.5/(ffit[ii+1]-ffit[ii]);
        }
       }
    
    //for(i=-(Nsten-1); i< 1; i++)
       for(i=-(smoothline->Nsten); i< 1; i++)
        {
         ii = iu+i;
         if(ii > -1 && ii < Nknot-1)
         {
             k = ii+(smoothline->Nsten-2);
             if(k > Nknot-2) k = Nknot-2;
             bb[ii] = 0.5*(slopes[k]+slopes[ii]);
         }
        }
    
    for(i=-2; i< 1; i++)
     {
      ii = iu+i;
      if(ii > -1 && ii < Nknot-2)
      {
          cc[ii+1] = 0.5*(slopes[ii+1]-slopes[ii]);
      }
     }
        
}

void sline(struct Smooth *smoothline, int istart, int iend, int Nknot, double *farray,  double *SM, double *ffit, double *slinep, double *bb, double *cc, double *dd)
{
    double f;
    int i, j, k, jx, jj, kk, ic;
    double x, y, z;
    double ax, bx;
    
    double *cx, *dx, *xv;
    
    cx = double_vector(smoothline->Nsten);
    dx = double_vector(smoothline->Nsten);
    xv = double_vector(smoothline->Nsten);
    
    ic = (smoothline->Nsten-1)/2;  // central point of stencil
    if(ic == 1) ic = 2;  // Nsten = 3 is a special case
    
    j = -ic+1;
   
    for (i = istart; i < iend; ++i)
    {
        f = farray[i];
        
          j -= 1;
          do
          {
           j++;
          } while(j < Nknot-smoothline->Nsten && f > 0.5*(ffit[j+ic]+ffit[j+ic-1]));
         
            
        jx = j;
        if(j < 0) jx = 0;
        if(j > Nknot-smoothline->Nsten) jx = Nknot-smoothline->Nsten;
        
        for(k=0; k< smoothline->Nsten-1; k++)
        {
            jj = jx+k;
            dx[k] = dd[jj];
            cx[k] = cc[jj];
            xv[k] = ffit[jj];
        }
        
        ax = slinep[jx];
        for(k=1; k< smoothline->Nsten-1; k++)
         {
            ax -= cx[k]*fabs((xv[k]-xv[0]));
        }
        bx = bb[jx];
                   
        y = ax + bx*(f-xv[0]);
        for(k=1; k< smoothline->Nsten-1; k++)
        {
            x = fabs(smoothline->asmooth*dx[k]*(f-xv[k]));
            if(x > 6.0)
            {
            y += cx[k]*fabs(f-xv[k]);
            }
            else
            {
             y += cx[k]/(smoothline->asmooth*dx[k])*ltc(smoothline,x);
            }
        }
        
        SM[i] = exp(y);
        
    }

     free_double_vector(cx);
     free_double_vector(dx);
     free_double_vector(xv);
    
    
}


void qscan(double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
       for (i = 0; i < N; ++i) dcopy[i] = data[i];

    // Tukey window
    tukey(dcopy, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(dcopy, 1, N);
    // whiten data
    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qtran.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}

//uses frequency domain data
void qscanf(double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;
    
    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];

    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    //TransformC(dcopy, freqs, tfD, tfDR, tfDI, 64, Tobs, N, Nf);
    
    
    out = fopen("Qsig.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}

//uses frequency domain data
void qscanshiftf(double **data, double **Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j, ii, k;
    double **tfDR, **tfDI, **tfD0, **tfD1;
    double fac, QQ;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
     // logarithmic frequency spacing
     subscale = 40;  // number of semi-tones per octave
     octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
     Nf = subscale*octaves+1;
     freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
     dx = log(2.0)/(double)(subscale);
     dlnf = dx;
     x = log(fmn);
     for(i=0; i< Nf; i++)
     {
         freqs[i] = exp(x);
         x += dx;
     }
     
    // printf("%d freqs %d times\n", Nf, N);
     
     tfDR = double_matrix(Nf,N);
     tfDI = double_matrix(Nf,N);
     tfD0 = double_matrix(Nf,N);
     tfD1 = double_matrix(Nf,N);
    
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    
    QQ = 16.0;
    
    for (i = 0; i < N/2; ++i) SX[i] = Sn[0][i]/fac;
    for (i = 0; i < N; ++i) dcopy[i] = data[0][i];
    whiten(dcopy, SX, N);
    TransformC(dcopy, freqs, tfD0, tfDR, tfDI, QQ, Tobs, N, Nf);
 
    for (i = 0; i < N/2; ++i) SX[i] = Sn[1][i]/fac;
    for (i = 0; i < N; ++i) dcopy[i] = data[1][i];
    whiten(dcopy, SX, N);
    TransformC(dcopy, freqs, tfD1, tfDR, tfDI, QQ, Tobs, N, Nf);
    
    k = (int)((9.877352e-04)/dt);
    
    out = fopen("Qdual.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            x = tfD0[j][i];
            ii = i-k;
            if(ii >= 0 && ii < N) x+= tfD1[j][i];
            fprintf(out,"%e %e %e\n", t, f, x);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    out = fopen("Q0.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            x = tfD0[j][i];
            fprintf(out,"%e %e %e\n", t, f, x);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    out = fopen("Q1.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            x = tfD1[j][i];
            fprintf(out,"%e %e %e\n", t, f, x);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD0,Nf);
    free_double_matrix(tfD1,Nf);
    
}


//uses frequency domain data
void qscanmaxf(double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j, k;
    double **tfDR, **tfDI, **tfD, **tfM;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];

    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    tfM = double_matrix(Nf,N);
    
    for(j = 0; j < Nf; j++)
    {
        for(i = 0; i < N; i++)
        {
            tfM[j][i] = 0.0;
        }
    }
    
    double Qx;
    
    Qx = 8.0;
    
    for(k = 0; k < 6; k++)
    {
        
    printf("%f\n", Qx);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qx, Tobs, N, Nf);
        
    Qx *= 2.0;
        
    for(j = 0; j < Nf; j++)
       {
           for(i = 0; i < N; i++)
           {
               if(tfD[j][i] > tfM[j][i]) tfM[j][i] = tfD[j][i];
           }
       }
    
    }
    
    
    out = fopen("Qmax.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfM[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    free_double_matrix(tfM,Nf);
    
}



void qscanres(double *data, double *signal, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn, fac;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window
    tukey(dcopy, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(dcopy, 1, N);
    // subtract the signal
    for (i = 0; i < N; ++i) dcopy[i] -= signal[i];
    // whiten data
    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qtranres.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}


void whiten(double *data, double *Sn, int N)
{
    double f, x, y, fix;
    int i;
    
    data[0] = 0.0;
    data[N/2] = 0.0;
    
    for(i=1; i< N/2; i++)
    {
        x = 1.0/sqrt(Sn[i]);
        data[i] *= x;
        data[N-i] *= x;
    }
    
}


void tukey(double *data, double alpha, int N)
{
  int i, imin, imax;
  double filter;
  
  imin = (int)(alpha*(double)(N-1)/2.0);
  imax = (int)((double)(N-1)*(1.0-alpha/2.0));
  
    int Nwin = N-imax;

 for(i=0; i< N; i++)
  {
    filter = 1.0;
    if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
    if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
    data[i] *= filter;
  }
  
}


void tukey_scale(double *s1, double *s2, double alpha, int N)
{
    int i, imin, imax;
    double x1, x2;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    int Nwin = N-imax;
    
    x1 = 0.0;
    x2 = 0.0;
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
        x1 += filter;
        x2 += filter*filter;
    }
    x1 /= (double)(N);
    x2 /= (double)(N);
    
    *s1 = x1;
    *s2 = sqrt(x2);
    
}



void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int imin, int imax, int N)
{
    int nb2, i, l, k, j;
    
    for (i = 0; i < N; i++)
    {
        corr[i]	= 0.0;
        corrf[i] = 0.0;
    }
    
    for (i=imin; i < imax; i++)
    {
        l=i;
        k=N-i;
        
        corr[l]	= (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]	= (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
        
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, N);
    
    
}


double fourier_nwip(double *a, double *b, double *Sn, int imin, int imax, int N)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=imin; i<imax; i++)
    {
        j = i;
        k = N-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}



void max_array_element(double *max, int *index, double *array, int n)
{
    int i;
    
    *max = array[0];
    *index = 0;
    
    for(i = 1; i <= n-1; i++)
    {
        if(array[i] > *max)
        {
            *max = array[i];
            *index = i;
        }
    }
}




void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n / 2;
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]    = (data1[l]*data2[l] + data1[k]*data2[k]);
        corr[k]    = (data1[k]*data2[l] - data1[l]*data2[k]);
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}


double det(double *A, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    int s;
    
    gsl_permutation *p = gsl_permutation_alloc(N);
    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[j*N+i]);
        }
    }
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    dx = 1.0;
    for (i = 0; i < N; i++) dx *= gsl_matrix_get (m, i, i);
    dx = fabs(dx);
    
    //returns the absolute value of the determinant.
    
    gsl_permutation_free(p);
    gsl_matrix_free(m);
    
    
    return dx;
    
    
}


double f_nwip(double *a, double *b, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-j;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product;
    }
    
    return(arg);
    
}

void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m)
{
    int j;
    double fix;
    
    // [0] t0 [1] f0 [2] Q [3] Amp [4] phi
    
    fix = sqrt((double)(n/2));
    
    #pragma omp parallel for
    for(j = 0; j < m; j++)
    {
        layerC(a, freqs[j], tf[j], tfR[j], tfI[j], Q, Tobs, fix, n);
    }
    
}


void layerC(double *a, double f, double *tf, double *tfR, double *tfI, double Q, double Tobs, double fix, int n)
{
    int i;
    double *AC, *AF;
    double *b;
    double *params;
    double bmag;
    
    params= double_vector(6);
    AC=double_vector(n);  AF=double_vector(n);
    b = double_vector(n);
    
    params[0] = 0.0;
    params[1] = f;
    params[2] = Q;
    params[3] = 1.0;
    params[4] = 0.0;
    
    SineGaussianC(b, params, Tobs, n);
    
    bmag = sqrt(f_nwip(b, b, n)/(double)n);
    
    bmag /= fix;
    
    phase_blind_time_shift(AC, AF, a, b, n);
    
    for(i = 0; i < n; i++)
    {
        tfR[i] = AC[i]/bmag;
        tfI[i] = AF[i]/bmag;
        tf[i] = tfR[i]*tfR[i]+tfI[i]*tfI[i];
    }
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(b);  free_double_vector(params);
    
}

void SineGaussianC(double *hs, double *sigpar, double Tobs, int N)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmx, fmn;//, fac;
    double phi, f;
    double tau;
    double re,im;
    double q, p, u;
    double A, B, C;
    double Am, Bm, Cm;
    double a, b, c;
    double dt, fac;
    double cosPhase_m, sinPhase_m, cosPhase_p, sinPhase_p;
    
    int i, imid, istart,istop,imin,imax,iend,even,odd;
    
    t0  = sigpar[0];
    f0  = sigpar[1];
    Q   = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmx = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmn = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    i = (int)(f0*Tobs);
    imin = (int)(fmn*Tobs);
    imax = (int)(fmx*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 1;
    if(imax > N/2) imax = N/2;
    
    hs[0] = 0.0;
    hs[N/2] = 0.0;
    
    for(i = 1; i < N/2; i++)
    {
        hs[i] = 0.0;
        hs[N-i] = 0.0;
    }
    
    dt = Tobs/(double)(N);
    fac = sqrt(sqrt(2.0)*PI*tau/dt);
    
    /* Use recursion relationship  */
    
    imid = (int)(f0*Tobs);
    
    p = PI*PI*tau*tau/(Tobs*Tobs);
    
    Bm = exp(-p*(((double)(imid)-f0*Tobs)*((double)(imid)-f0*Tobs)));
    Cm = 1.0;
    
    b = exp(-p*(1.0+2.0*((double)(imid)-f0*Tobs)));
    c = exp(-2.0*p);
    
    // start in the middle and work outwards
    
    B = Bm;
    C = Cm;

    
    for(i = imid; i < imax; i++)
    {
 
        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        //  printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        
    }
    
    // reset to midpoint
    
    B = Bm;
    C = Cm;
    
    b = exp(p*(-1.0+2.0*((double)(imid)-f0*Tobs)));
    // c unchanged
    
    for(i = imid; i > imin; i--)
    {

        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        // printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        

    }
    
    
}

void Inverse(double **M, double **IM, int d)
{
    int i, j;
    int s;
    double x, maxc;
    
    gsl_matrix *m = gsl_matrix_alloc (d, d);
    
    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            gsl_matrix_set(m, i, j, M[i][j]);
        }
    }
    
    gsl_permutation *p = gsl_permutation_alloc(d);
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(d, d);
    
    gsl_linalg_LU_invert(m, p, inv);
    
    gsl_permutation_free(p);
    
    
    for (i = 0; i < d; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            IM[i][j] = gsl_matrix_get(inv, j, i);
        }
    }
    
    
    gsl_matrix_free (inv);
    gsl_matrix_free (m);
    
    return;
    
}


void makespec(int Nspline, int Nlines, double *ffit, double *Sspline, double *linef, double *lineh,  double *lineQ, double *linew, double *deltafmax, double *SN, double Tobs,  int N)
{

    int i, j;
    double f, x, y, z;
    double sm;
    double f2, ff2, df2, df4;
    
    // Allocate spline
    gsl_spline   *cspline = gsl_spline_alloc(gsl_interp_cspline, Nspline);
    gsl_interp_accel *acc    = gsl_interp_accel_alloc();
    
    /* compute spline */
    gsl_spline_init(cspline,ffit,Sspline,Nspline);
    
    SN[0] = exp(Sspline[0]);
    SN[N/2-1] = exp(Sspline[Nspline-1]);
    
    for (i = 1; i < N/2-1; ++i)
       {
           
           f = (double)(i)/Tobs;
           
           sm = exp(gsl_spline_eval(cspline,f,acc));
           
           y = 0.0;
           for (j = 0; j < Nlines; ++j)
           {
              
               x = fabs(f - linef[j]);
               
               if(x < linew[j])
               {
                   z = 1.0;
                   if(x > deltafmax[j]) z = exp(-(x-deltafmax[j])/deltafmax[j]);
                   f2 = linef[j]*linef[j];
                   ff2 = f*f;
                   df2 = lineQ[j]*(1.0-ff2/f2);
                   df4 = df2*df2;
                   y += z*lineh[j]/(ff2/f2+df4);
               }
           }
          
           SN[i] = sm*(1.0+y);
           
       }
    
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc);
    
}



void specest(double *data, double *Hf, int N, int Ns, double dt, double fmx, double *SN, double *SM, double *PS)
{
    int i, j, k, M, Nf, Nstep,  ii, m;
    int jj, kk, Nlines;
    int oflag, flag;
    int imin, imax;
    double SNR, max;
    double junk, Tobs, fix, f, t, t0, df, x, y, z, dx;
    double fmn, dfx, Q, fny, delt, scale, dlnf;
    double *freqs, *ref;
    double *Draw;
    double *D, *times;
    double *Sn;
    double *specD, *sspecD;
    double *sdata;
    double *intime, *sqf;
    double *pspline, *dspline;
    double sigmean, sigmedian;
    int subscale, octaves;
    int mmax;
    double SNRsq, SNRold, pH, pL, pmax;
    double SNRH, SNRL, pw, alpha;
    double t_rise, s1, s2, ascale, fac, Dfmax;
    double *linew;
    int pflag;
    
    int modelprint;

    char filename[1024];
    char command[1024];
    
    FILE *out;
    
    Q = Qs;  // Q of transform
    Dfmax = 8.0;  // width of smoothing window in Hz
    Tobs = (double)(N)*dt;  // duration
    
    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    
    
    df = 1.0/Tobs;  // frequency resolution
    fny = 1.0/(2.0*dt);  // Nyquist
    
    // Choose the range of the spectrogram. fmin, fmax must be a power of 2
    fmn = 8.0;
    if(fmx > fny) fmx = fny;
    
    D = (double*)malloc(sizeof(double)* (N));
    Draw = (double*)malloc(sizeof(double)* (N));
    
    // Copy data over
    for(i=0; i< N; i++)
    {
     Draw[i] = data[i];
     D[i] = data[i];
    }
    
    // Normalization factor
    tukey_scale(&s1, &s2, alpha, N);
    
    printf("%f %f\n", s1, s2);
    
    // Time series data is corrupted by the Tukey window at either end
    // imin and imax define the "safe" region to use
    imin = (int)(2.0*t_rise/dt);
    imax = N-imin;
    
    
    // Prepare to make spectogram
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    sqf = (double*)malloc(sizeof(double)* (Nf));
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        sqf[i] = sqrt(freqs[i]);
        x += dx;
    }
    
    
    
    printf("frequency layers = %d\n", Nf);
    
    scale = Getscale(freqs, Q, Tobs, fmx, N, Nf);
    
    sspecD = (double*)malloc(sizeof(double)*(N/2));
    specD = (double*)malloc(sizeof(double)*(N/2));
    Sn = (double*)malloc(sizeof(double)*(N/2));
    

    SNRold = 0.0;
    clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Q, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 0);
    
    
    // if big glitches are present we need to rinse and repeat
    i = 0;
    while(i < 10 && (SNR-SNRold) > 10.0)
    {
        SNRold = SNR;
        clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Q, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 0);
        i++;
    }
    
    
    // pass back the cleaned data. This cleaned data is not passed forward
    // from the spectral estimation code. Only used inside the MCMC
    for(i=0; i< N; i++)
    {
     data[i] = D[i];
    }
    
    // re-compute the power spectrum using the cleaned data
    tukey(D, alpha, N);
    
    gsl_fft_real_radix2_transform(D, 1, N);
    
    // Form spectral model for whitening data (lines plus a smooth component)
    spectrum(D, Sn, specD, sspecD, df, N);
    
    fac = Tobs/((double)(N)*(double)(N));
    
    for (i = 0; i < Ns; ++i) SM[i] = sspecD[i]*fac;
    for (i = 0; i < Ns; ++i) SN[i] = specD[i]*fac;
    for (i = 0; i < Ns; ++i) PS[i] = Sn[i]*fac;
    
    // print the cleaned data
    for(i=0; i< N; i++)
    {
     D[i] = data[i];
    }
    
    // save Q-scan data to file
    if(verbose == 1)
    {
    clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Qprint, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 1);
    }
    
    free(D);
    free(Draw);
    free(sspecD);
    free(specD);
    free(Sn);
    free(freqs);
    free(sqf);

    
}


void clean(double *D, double *Draw, double *Hf, double *sqf, double *freqs, double *Sn, double *specD, double *sspecD, double df, double Q, double Tobs, double scale, double alpha, int Nf, int N, int imin, int imax, double *SNR, int pflag)
{
    
    int i, j, k;
    int flag;
    int ii, jj;
    double x, y, dt;
    double S;
    double fac;
    double t, f, fmn;
    
    // allocate some arrays
    double **tfDR, **tfDI;
    double **tfD;
    double **live;
    double **live2;
    
    double *Dtemp, *Drf;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *out;
    
    live= double_matrix(Nf,N);
    live2= double_matrix(Nf,N);
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    dt = Tobs/(double)(N);
    
    
    Dtemp = (double*)malloc(sizeof(double)*(N));
    
    for (i = 0; i < N; i++) Dtemp[i] = Draw[i];
    
    // D holds the previously cleaned time-domain data
    // Draw holds the raw time-domain data
    
    // D is used to compute the spectrum. A copy of Draw, Dtemp, is then whitened using this spectrum and glitches are then identified
    // The glitches are then re-colored, subtracted from Draw to become the new D
 
 
    // Tukey window
    tukey(D, alpha, N);
    tukey(Dtemp, alpha, N);
   
    
    // FFT
    gsl_fft_real_radix2_transform(D, 1, N);
    gsl_fft_real_radix2_transform(Dtemp, 1, N);
    
    // remove the CBC signal if provided (in the spec code this is zero)
    for(i = 0; i < N; i++)
    {
        D[i] -= Hf[i];
        Dtemp[i] -= Hf[i];
    }
    
    // Form spectral model for whitening data (lines plus a smooth component)
       start = clock();
       spectrum(D, Sn, specD, sspecD, df, N);
       end = clock();
       cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
       printf("spectrum took %f seconds\n", cpu_time_used);
   
    
    // whiten data
    whiten(Dtemp, specD, N);
    
    // Wavelet transform
    start = clock();
    TransformC(Dtemp, freqs, tfD, tfDR, tfDI, Q, Tobs, N, Nf);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("wavelet transform took %f seconds\n", cpu_time_used);
    
    
    if(pflag == 1)
    {
        out = fopen("Qtransform.dat","w");
        for(j = 0; j < Nf; j++)
        {
            f = freqs[j];
            
            for(i = 0; i < N; i++)
            {
                t = (double)(i)*dt-Tobs+2.0;  // trigger time is set two seconds from end, set zero there
                fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
            }
            
            fprintf(out,"\n");
        }
        fclose(out);
        
    }
    
    k = 0;
    //  apply threshold
    for(j = 0; j < Nf; j++)
    {
        for (i = 0; i < N; i++)
        {
            live[j][i] = -1.0;
            if(tfD[j][i] > sthresh) live[j][i] = 1.0;
            if(tfD[j][i] > sthresh) k++;
            live2[j][i] = live[j][i];
        }
    }
   
    
    // dig deeper to extract clustered power
    for(j = 1; j < Nf-1; j++)
    {
        for (i = 1; i < N-1; i++)
        {
            
            flag = 0;
            for(jj = -1; jj <= 1; jj++)
            {
                for(ii = -1; ii <= 1; ii++)
                {
                    if(live[j+jj][i+ii] > 0.0) flag = 1;
                }
            }
            if(flag == 1 && tfD[j][i] > warm) live2[j][i] = 1.0;
        }
    }
    
    
    // build the excess power model
    for (i = 0; i < N; i++)
    {
        Dtemp[i] = 0.0;
    }
    
    k = 0;
    for(j = 0; j < Nf; j++)
    {
        for (i = imin; i < imax; i++)
        {
            if(live2[j][i] > 0.0) Dtemp[i] += scale*sqf[j]*tfDR[j][i];
        }
    }
    
    if(pflag == 1)
       {
        out = fopen("wglitch.dat", "w");
        for(i=0; i< N; i++)
        {
            fprintf(out,"%e %e\n", (double)(i)*dt-Tobs+2.0, Dtemp[i]);
         }
       fclose(out);
       }
    
    // Compute the excess power (relative to the current spectral model
    S = 0.0;
    for (i = imin; i < imax; ++i) S += Dtemp[i]*Dtemp[i];
    S = sqrt(S);
    
   
    printf("Excess SNR at Q %f = %f\n", Q, S);
    
    
    *SNR = S;
    

    //Unwhiten and subtract the excess power so we can compute a better spectral estimate
    // Back to frequency domain
    
    gsl_fft_real_radix2_transform(Dtemp, 1, N);
    
    
    fmn = freqs[0];
    
    // only use smooth spectrum in the un-whitening
    Dtemp[0] = 0.0;
    for(i=1; i< N/2; i++)
    {
        f = (double)(i)/Tobs;
        //y = 1.0;
        //if(f < fmn) y = 0.5*(1.0+tanh(8.0*(f-0.5*fmn)/fmn));
        //x = y*sqrt(sspecD[i]);
        x = sqrt(sspecD[i]);
        Dtemp[i] *= x;
        Dtemp[N-i] *= x;
    }
    Dtemp[N/2] = 0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(Dtemp, 1, N);
    
    
    x = sqrt((double)(2*N));
    
    // note that this D[i] still contains the signal since it
    // is formed from the raw data minus the glitch model
    for(i=0; i< N; i++)
    {
        D[i] = Draw[i]-Dtemp[i]/x;
    }
    
    /*
    out = fopen("glitch.dat", "w");
    for(i=0; i< N; i++)
    {
        fprintf(out,"%d %e %e\n", i, Dtemp[i]/x, Draw[i]);
    }
    fclose(out);
    */
    
    free(Dtemp);
    free_double_matrix(live,Nf);
    free_double_matrix(live2,Nf);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
    
    return;
    
    
}

void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, int N)
{
    double Df, Dfmax, x, y;
    double Df1, Df2;
    int mw, k, i, j;
    int mm, kk;
    int end1, end2, end3;
    double med;
    double *chunk;
    
    // log(2) is median/2 of chi-squared with 2 dof
    
    
    for(i=1; i< N/2; i++) S[i] = 2.0*(data[i]*data[i]+data[N-i]*data[N-i]);
    S[0] = S[1];
    
    Dfmax = 16.0; // is the  width of smoothing window in Hz
    
    // Smaller windows used initially where the spectrum is steep
    Df2 = Dfmax/2.0;
    Df1 = Dfmax/4.0;
    
    // defines the ends of the segments where smaller windows are used
    end1 = (int)(32.0/df);
    end2 = 2*end1;
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    //printf("numer of bins in smoothing window %d\n", mw);
    k = (mw+1)/2;
    chunk = double_vector(mw);
    
    end3 = N/2-k;  // end of final chunk
    
    // Fill the array so the ends are not empty - just to be safe
    for(i=0;i< N/2;i++)
    {
        Sn[i] = S[i];
        Smooth[i] = S[i];
    }
    

    mw = (int)(Df1/df)+1;  // size of median window
    k = (mw+1)/2;
    
    for(i=4; i< k; i++)
    {
        mm = i/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        gsl_sort(chunk, 1, mm);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mm)/LN2;
    
        Smooth[i] = Sn[i];
        
    }
    
    
    i = k;
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
        
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end1);
    
    
    
    mw = (int)(Df2/df)+1;  // size of median window
    k = (mw+1)/2;
    
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
        
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end2);
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    k = (mw+1)/2;
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
    
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end3);
    
    
    for(i=end3; i< N/2-4; i++)
    {
        mm = (N/2-i)/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        gsl_sort(chunk, 1, mm);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mm)/LN2;
        
        Smooth[i] = Sn[i];
        
    }
    
    
    free_double_vector(chunk);
    
    
    
    
    // zap the lines.
    for(i=1;i< N/2;i++)
    {
        x = S[i]/Sn[i];
        if(x > linemul)
        {
            Sn[i] = S[i];
        }
    }
    
    
    
}



void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2)
{
    /* Butterworth bandpass filter
     n = filter order 4,8,12,...
     s = sampling frequency
     f1 = upper half power frequency
     f2 = lower half power frequency  */
    
    if(n % 4){ printf("Order must be 4,8,12,16,...\n"); return;}
    
    int i, j;
    double a = cos(PI*(f1+f2)/s)/cos(PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    
    n = n/4;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *d3 = (double *)malloc(n*sizeof(double));
    double *d4 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)malloc(n*sizeof(double));
    double *w1 = (double *)malloc(n*sizeof(double));
    double *w2 = (double *)malloc(n*sizeof(double));
    double *w3 = (double *)malloc(n*sizeof(double));
    double *w4 = (double *)malloc(n*sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*(double)i+1.0)/(4.0*(double)n));
        s = b2 + 2.0*b*r + 1.0;
        A[i] = b2/s;
        d1[i] = 4.0*a*(1.0+b*r)/s;
        d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        d3[i] = 4.0*a*(1.0-b*r)/s;
        d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
        w3[i] = 0.0;
        w4[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + x;
            x = A[i]*(w0[i] - 2.0*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
    
    return;
}




double Getscale(double *freqs, double Q, double Tobs, double fmx, int n, int m)
{
    double *data, *intime, *ref, **tfR, **tfI, **tf;
    double f, t0, delt, t, x, fix, dt;
    double scale, sqf;
    int i, j;
    
    FILE *out;
    
    data = double_vector(n);
    ref = double_vector(n);
    intime = double_vector(n);
    
    tf = double_matrix(m,n);
    tfR = double_matrix(m,n);
    tfI = double_matrix(m,n);
    
    f = fmx/4.0;
    t0 = Tobs/2.0;
    delt = Tobs/8.0;
    dt = Tobs/(double)(n);
    
    //out = fopen("packet.dat","w");
    for(i=0; i< n; i++)
    {
        t = (double)(i)*dt;
        x = (t-t0)/delt;
        x = x*x/2.0;
        data[i] = cos(TPI*t*f)*exp(-x);
        ref[i] = data[i];
        //fprintf(out,"%e %e\n", t, data[i]);
    }
    // fclose(out);
    
    gsl_fft_real_radix2_transform(data, 1, n);

    TransformC(data, freqs, tf, tfR, tfI, Q, Tobs, n, m);
    
    for(i = 0; i < n; i++) intime[i] = 0.0;
    
    for(j = 0; j < m; j++)
    {
        
        f = freqs[j];
        
        
         sqf = sqrt(f);
        
        for(i = 0; i < n; i++)
        {
            intime[i] += sqf*tfR[j][i];
        }
        
    }
    
    x = 0.0;
    j = 0;
    // out = fopen("testtime.dat","w");
    for(i=0; i< n; i++)
    {
        // fprintf(out,"%e %e %e\n",times[i], intime[i], ref[i]);
        
        if(fabs(ref[i]) > 0.01)
        {
            j++;
            x += intime[i]/ref[i];
        }
    }
    //fclose(out);
    
    x /= sqrt((double)(2*n));
    
    scale = (double)j/x;
    
    // printf("scaling = %e %e\n", x/(double)j, (double)j/x);
    
    free_double_vector(data);
    free_double_vector(ref);
    free_double_vector(intime);
    free_double_matrix(tf,m);
    free_double_matrix(tfR,m);
    free_double_matrix(tfI,m);
    
    return scale;
    
}


void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase)
{
    /* Update re and im for the next iteration. */
    double cosphi = *cosPhase;
    double sinphi = *sinPhase;
    double x, y;
    
    x = (cosphi*dre + sinphi*dim);
    y = (sinphi*dre - cosphi*dim);
    
    double newRe = cosphi - x;
    double newIm = sinphi - y;
    
    *cosPhase = newRe;
    *sinPhase = newIm;
    
}

void SineGaussianF(double *hs, double *sigpar, double Tobs, int N)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmx, fmn;//, fac;
    double phi, f;
    double tau;
    double re,im;
    double q, p, u;
    double A, B, C;
    double Am, Bm, Cm;
    double a, b, c;
    double cosPhase_m, sinPhase_m, cosPhase_p, sinPhase_p;
    
    int i, imid, istart,istop,imin,imax,iend,even,odd;
    
    t0  = sigpar[0];
    f0  = sigpar[1];
    Q   = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmx = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmn = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    i = (int)(f0*Tobs);
    imin = (int)(fmn*Tobs);
    imax = (int)(fmx*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 1;
    if(imax > N/2) imax = N/2;
    
    hs[0] = 0.0;
    hs[N/2] = 0.0;
    
    for(i = 1; i < N/2; i++)
    {
        hs[i] = 0.0;
        hs[N-i] = 0.0;
    }
    
    /* Use recursion relationship  */
    
    //incremental values of exp(iPhase)
    double dim = sin(TPI*t0/Tobs);
    double dre = sin(0.5*(TPI*t0/Tobs));
    dre = 2.0*dre*dre;
    
    double amplitude = 0.5*(Amp)*RTPI*tau;
    double pi2tau2   = PI*PI*tau*tau;
    double Q2        = Q*Q/f0;
    
    
    imid = (int)(f0*Tobs);
    
    
    q = Q*Q/(f0*Tobs);
    p = PI*PI*tau*tau/(Tobs*Tobs);
    u = PI*PI*tau*tau/Tobs*(2.0*f0*Tobs-1.0);
    
    Am = exp(-q*(double)(imid));
    Bm = exp(-p*(((double)(imid)-f0*Tobs)*((double)(imid)-f0*Tobs)));
    Cm = 1.0;
    
    a = exp(-q);
    b = exp(-p*(1.0+2.0*((double)(imid)-f0*Tobs)));
    c = exp(-2.0*p);
    
    // sine and cosine of phase at reference frequency
    f = (double)(imid)/Tobs;
    double phase = TPI*f*t0;
    double cosPhase_m0  = cos(phase-phi);
    double sinPhase_m0  = sin(phase-phi);
    double cosPhase_p0  = cos(phase+phi);
    double sinPhase_p0  = sin(phase+phi);
    
    // start in the middle and work outwards
    
    A = Am;
    B = Bm;
    C = Cm;
    
    cosPhase_m  = cosPhase_m0;
    sinPhase_m  = sinPhase_m0;
    cosPhase_p  = cosPhase_p0;
    sinPhase_p  = sinPhase_p0;
    
    for(i = imid; i < imax; i++)
    {
        even = 2*i;
        odd = even + 1;
        f = (double)(i)/Tobs;
        
        sf = amplitude*B;
        sx = A;
        re = sf*(cosPhase_m+sx*cosPhase_p);
        im = -sf*(sinPhase_m+sx*sinPhase_p);
        
        //  printf("%e %e\n", exp(-Q2*f)-A, exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = re;
        hs[N-i] = im;
        
        A *= a;
        B *= (C*b);
        C *= c;
        
        /* Now update re and im for the next iteration. */
        recursive_phase_evolution(dre, dim, &cosPhase_m, &sinPhase_m);
        recursive_phase_evolution(dre, dim, &cosPhase_p, &sinPhase_p);
    }
    
    // reset to midpoint
    
    A = Am;
    B = Bm;
    C = Cm;
    
    cosPhase_m  = cosPhase_m0;
    sinPhase_m  = sinPhase_m0;
    cosPhase_p  = cosPhase_p0;
    sinPhase_p  = sinPhase_p0;
    
    a = 1.0/a;
    b = exp(p*(-1.0+2.0*((double)(imid)-f0*Tobs)));
    // c unchanged
    
    // interate backwards in phase
    dim *= -1.0;
    
    for(i = imid; i > imin; i--)
    {
        even = 2*i;
        odd = even + 1;
        f = (double)(i)/Tobs;
        
        sf = amplitude*B;
        sx = A;
        re = sf*(cosPhase_m+sx*cosPhase_p);
        im = -sf*(sinPhase_m+sx*sinPhase_p);
        
        // printf("%e %e\n", exp(-Q2*f)-A, exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = re;
        hs[N-i] = im;
        
        A *= a;
        B *= (C*b);
        C *= c;
        
        /* Now update re and im for the next iteration. */
        recursive_phase_evolution(dre, dim, &cosPhase_m, &sinPhase_m);
        recursive_phase_evolution(dre, dim, &cosPhase_p, &sinPhase_p);
    }
    
    
}



void shift(double *a, double R, double delt, double pshift, int n, double Tobs)
{
    int i, j, k;
    double  f;
    double ReA, ImA;
    double cx, sx;
    
    for(i=0; i<n/2; i++)
    {
        f = (double)(i)/Tobs;
        j = i; // real
        k = n-i; // imaginary
        cx = cos(-pshift+TPI*delt*f);
        sx = sin(-pshift+TPI*delt*f);
        ReA = a[j]*cx-a[k]*sx;
        ImA = a[k]*cx+a[j]*sx;
        a[j] = R*ReA;
        a[k] = R*ImA;
    }
    
    return;
    
}


double fourier_nwip_shift(double *a, double *b, double delt, double pshift, int n, double Tobs, int imin, int imax)
{
    int i, j, k;
    double arg, product, f;
    double ReA, ReB, ImA, ImB;
    double cx, sx;
    
    // Does f_nwip with a given time and phase shift
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        if(i > imin && i < imax)
        {
            f = (double)(i)/Tobs;
            j = i; // real
            k = n-i; // imaginary
            ReA = a[j];
            ImA = a[k];
            cx = cos(-pshift+TPI*delt*f);
            sx = sin(-pshift+TPI*delt*f);
            ReB = b[j]*cx-b[k]*sx;
            ImB = b[k]*cx+b[j]*sx;
            product = ReA*ReB + ImA*ImB;
            arg += product;
        }
    }
    
    return(arg);
    
}






int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}

int **int_matrix(int N, int M)
{
    int i;
    int **m = malloc( (N+1) * sizeof(int *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(int));
    }
    
    return m;
}

void free_int_matrix(int **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_int_vector(m[i]);
    free(m);
}

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    int i;
    double **m = malloc( (N+1) * sizeof(double *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(double));
    }
    
    return m;
}

void free_double_matrix(double **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_double_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
        }
    }
    
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
    
    free(t);
}

double ****double_quad(int N, int M, int L, int K)
{
    int i,j,k;
    
    double ****t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
            for(k=0; k<L+1; k++)
            {
                       t[i][j][k] = malloc( (K+1) * sizeof(double));
            }
        }
    }
    
    return t;
}

void free_double_quad(double ****t, int N, int M, int L)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_tensor(t[i],M,L);
    
    free(t);
}







