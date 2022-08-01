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
#include <ctype.h>


#define NR_END 1
#define FREE_ARG char*

/* gcc -o getspec_gwosc getspec_gwosc.c -lm */


int main(int argc, char *argv[])
{
  int i, j, k, n, m, M, N;
  double fmax;
  double ttrig, tstart, starttime, endtime;
  double Tobs, tbuf;
  double t, t0, dt, x;
	
  char filename[1024];
  char command[1024];
  char channel[1024];
  char obs[16];
    char line1[100];
    char line2[100];
    char line3[100];

  FILE *in;
  FILE *ifp;
  FILE *out;
    
    if(argc!=7)
    {
        printf("./gwoscdump filename Tobs trig_time time_to_end fmax observatory\n");
        printf("\n Tobs must be a multiple of 2s\n");
        printf("\n GPS trigger time gets rounded to nearest integer\n");
        printf("\n time_to_end is time between trigger and end of buffer in integer seconds (max is Tobs/2)\n");
        printf("\n maximum frequency for the analysis\n");
        printf("\n Observatory is H, L or V. Geo and Kagra not yet supported\n");
        return 1;
    }
    
    // Designed to use data from
    // https://www.gw-openscience.org/eventapi/html/allevents/

    // The code is happy reading in the 4 kHz or 16 kHz data in .txt format. 32 seconds is usually enough
    
    in = fopen(argv[1],"r");
    if (in == NULL){
        printf("Could not open file %s", argv[1]);
        return 2;
    }
    
    Tobs = atof(argv[2]);
    if((int)(Tobs)%2 != 0)
       {
           printf("Observation time must be a multiple of 2s\n");
           return 1;
       }
    
    ttrig = atof(argv[3]);
    
    tbuf = floor(atof(argv[4]));
    if(tbuf > Tobs/2.0) tbuf = Tobs/2.0;
    
    fmax = atof(argv[5]);
    
    sprintf(obs,"%s",argv[6]);
    k = strcmp(obs,"H");
    if(k==0) m = 0;
    k = strcmp(obs,"L");
    if(k==0) m = 1;
    k = strcmp(obs,"V");
    if(k==0) m = 2;
   
    tstart = ttrig + tbuf - Tobs; //Puts the trigger tbuf from the end of the data segment

    starttime = tstart;
    endtime = tstart + Tobs;
    
    // Strip off text header
    
    fgets(line1, 100, in);
    fgets(line2, 100, in);
    fgets(line3, 100, in);
    
    /* printf("%s\n", line1);
    printf("%s\n", line2);
    printf("%s\n", line3); */
    
    long num[4];
    
    
    char *ptr = line2;
    
    i = 0;
    while (*ptr)
    {
        if  (isdigit (*ptr) ){
            long  val = strtol (ptr, &ptr, 10);
            num[i] = val;
            i++;
        }  else {
            ptr++;
        }
    }
    
    char *pt = line3;
    
    while (*pt)
    {
        if  (isdigit (*pt) ){
            long  val = strtol (pt, &pt, 10);
            num[i] = val;
            i++;
        }  else {
            pt++;
        }
    }
    
    printf("cadence %ld GPS start %ld duration %ld\n", num[0], num[1], num[2]);


    N = (int)num[0]*(int)num[2];
    t0 = (double)num[1];
    
    dt = (double)(num[2])/(double)(N);
    
    if(tstart < t0 || endtime > t0+(double)(num[2]))
    {
        printf("Requested time range outside of supplied GWOSC data file\n");
    }
    
    sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
        

        out = fopen(command, "w");
        for(i=0; i< N; i++)
        {
            fscanf(in,"%lf", &x);
            t = t0 + (double)(i)*dt;
            if(t >= tstart && t < endtime)
            {
                fprintf(out,"%.16e %.16e\n",  t, x);
            }
        }
        fclose(out);
        fclose(in);
    
       sprintf(command, "./SpecFit %d %d %d %d", (int)(Tobs), (int)ttrig, m, (int)fmax);
       printf("%s\n", command);
       system(command);
        
       sprintf(command, "cp PSD.dat PSD_%d_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m, (int)fmax);
       printf("%s\n", command);
       system(command);
    
        sprintf(command, "cp white.dat white_%d_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m, (int)fmax);
        printf("%s\n", command);
        system(command);
    
        sprintf(command, "./Gauss_Test_Freq white.dat 8");
        printf("%s\n", command);
        system(command);
    
        sprintf(command, "gnuplot ADcheck.gnu");
        printf("%s\n", command);
        system(command);
    
    
    return 0;

}





