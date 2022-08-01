#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>


#define NR_END 1
#define FREE_ARG char*

static REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS
                                    start, REAL8 length);

static void output_frame(REAL8TimeSeries *timeData, CHAR *frameType, CHAR
                         *channel, CHAR *ifo);

int main(int argc, char *argv[])
{
  int i, j, k, M, N, Nf, Nstep, Nclean, ii, m, rs, tsi, tti;
    int jj, kk;
  int linlog;
  int imin, imax;
  double SNR;
  double junk, Tobs, fix, f, t, t0, dt, dtm, df, x, y, dx;
  double dfx, Q, fny, delt, scale, dlnf;
  double Hmax, Lmax;
  double pshift;
  double *freqs, *data, *ref;
  double *inp, *oup, *slice;
  double *tempH1;
  double *tempL1;
  double *H1dat, *L1dat;
  double *Hclean, *Lclean;
  double *Hraw, *Lraw, *traw;
  double *H, *L, *times;
  double *SH, *SL;
  double **geo;
  double *specH1, *specL1, *sspecH1, *sspecL1;
  double *sdata;
  double **tfmap, **tfavH1, **tfavL1;
  double *intime, *sqf;
  double sthresh, nthresh;
  double sigmean, sigmedian;
  int subscale, octaves;
    int mmax, flag;
    double SNRsq, pH, pL, pmax;
    double SNRH, SNRL, pw, alpha;
    int specprint, printglitch;
   double t_rise, s1, s2, ascale, warm;
  double ttrig, tbuf, tstart, tstart_clean, Tclean, starttime, endtime;
    int cflag, size;
    int fmax;
	
  char filename[1024];
  char command[1024];
  char L1name[1024];
  char H1name[1024];
    
    char channel[1024];
    char obs[16];
    char version[16];
    char ftype[64];
    
  int n;


  FILE *in;
  FILE *ifp;
  FILE *out;

    
    dt = 1.0/16384.0;
    
    if(argc!=8)
    {
        printf("\n Usage ./getspec Tobs trig_time time_to_end fmax observatory frame_version cleaned?\n");
        printf("\n Tobs must be a multiple of 2s\n");
        printf("\n GPS trigger time gets rounded to nearest integer\n");
        printf("\n time_to_end is time between trigger and end of buffer in integer seconds (max is Tobs/2)\n");
        printf("\n maximum frequency for the analysis)\n");
        printf("\n Observatory is H, L or V. Geo and Kagra not yet supported\n");
        printf("\n frame_version is typically C00, C01, C02. Also can be low latency type\n");
        printf("\n cleaned? is a flag that is either 0 or 1. Set to 1 if you want cleaned data\n");
        printf("\n You Can't Always Get What You Want. If you ask for say C02 cleaned, it might not exist\n");
        return 1;
    }
    
    Tobs = atof(argv[1]);
    ttrig = floor(atof(argv[2]));
    tbuf = floor(atof(argv[3]));
    if(tbuf > Tobs/2.0) tbuf = Tobs/2.0;
    fmax = atoi(argv[4]);
    sprintf(obs,"%s",argv[5]);
    sprintf(version,"%s",argv[6]);
    cflag = atoi(argv[7]);
    tstart = ttrig + tbuf - Tobs; //Puts the trigger tbuf from the end of the data segment
    
    if((int)(Tobs)%2 != 0)
    {
        printf("Observation time must be a multiple of 2s\n");
        return 1;
    }
    

    if(cflag == 0)  // not using cleaned frames
    {
        
    sprintf(ftype,"%s1_HOFT_%s",argv[5],argv[6]);
    k = strcmp(version,"C00");
    if(k==0)
    {
    sprintf(channel,"%s1:DCS-CALIB_STRAIN",argv[5]);
    }
    else
    {
     sprintf(channel,"%s1:DCS-CALIB_STRAIN_%s",argv[5],argv[6]);
    }
        
    }
    else  // using cleaned frames (there are no C00 cleaned frames)
    {
        sprintf(ftype,"%s1_CLEANED_HOFT_%s",argv[4],argv[5]);
        sprintf(channel,"%s1:GDS-CALIB_STRAIN_CLEAN",argv[4]);
    }
    
    // Virgo is a different beast
    k = strcmp(obs,"V");
    if(k==0)
     {
         if(cflag == 0)  // not using cleaned frames
         {
            sprintf(ftype,"V1Online");
            sprintf(channel,"V1:Hrec_hoft_16384Hz");
         }
         else
         {
          sprintf(ftype,"V1O2Repro2A");
          sprintf(channel,"V1:Hrec_hoft_V1O2Repro2A_16384Hz");
         }
     }
    
    printf("%s %s\n", ftype, channel);
    
    starttime = tstart;
    endtime = tstart + Tobs;
    
    N = (int)(Tobs/dt);
    

 
        sprintf(command, "gw_data_find --observatory %s --type %s -s %d -e %d --lal-cache -u file > f.cache", obs, ftype, (int)starttime, (int)endtime);
        printf("%s\n", command);
        system(command);
    
    // now check to see if anything was found
    in = fopen("f.cache","r");
    fseek (in, 0, SEEK_END);
    size = ftell(in);
    fclose(in);
    if (0 == size)
    {
        printf("Requested data not found. Check GPS time, frame and channel request\n");
        return -1;
    }
   
    REAL8TimeSeries* timeData=NULL;
    char * cachefile="f.cache";
        
    LIGOTimeGPS epoch;
        
    // Set time & retrieve data by reading frame cache
    XLALGPSSetREAL8(&epoch, tstart);
    
    printf("reading in the LIGO data\n");
    timeData=readTseries(cachefile,channel,epoch,Tobs);
    
    /*
    for(i=0; i< N; i++)
    {
        t = tstart + (double)(i)*dt;
        printf("%.16e %.16e\n",  t, timeData->data->data[i]);
    }
    */
    
        dt= timeData->deltaT;
        df = 1.0/Tobs;  // frequency resolution
        N = timeData->data->length;
    
        k = strcmp(obs,"H");
        if(k==0) m = 0;
        k = strcmp(obs,"L");
        if(k==0) m = 1;
        k = strcmp(obs,"V");
        if(k==0) m = 2;
    
        sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    
      
        out = fopen(command, "w");
        for(i=0; i< N; i++)
        {
            t = tstart + (double)(i)*dt;
            fprintf(out,"%.16e %.16e\n",  t, timeData->data->data[i]);
        }
        fclose(out);
    
    
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


/* ********************************************************************************** */
/*																					  */
/*                                 Frame I/O                                          */
/*																					  */
/* ********************************************************************************** */

static REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
    LALStatus status;
    memset(&status,0,sizeof(status));
    LALCache *cache = NULL;
    LALFrStream *stream = NULL;
    REAL8TimeSeries *out = NULL;
    
    cache  = XLALCacheImport( cachefile );
    int err;
    err = *XLALGetErrnoPtr();
    if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file \"%s\",\n       XLALError: \"%s\".\n",cachefile, XLALErrorString(err)); exit(-1);}
    stream = XLALFrStreamCacheOpen( cache );
    if(stream==NULL) {fprintf(stderr,"ERROR: Unable to open stream from frame cache file\n"); exit(-1);}
    out = XLALFrStreamInputREAL8TimeSeries( stream, channel, &start, length , 0 );
    if(out==NULL) fprintf(stderr,"ERROR: unable to read channel %s from %s at time %i\nCheck the specified data duration is not too long\n",channel,cachefile,start.gpsSeconds);
    XLALDestroyCache(cache);
    LALFrClose(&status,&stream);
    return out;
}


static void output_frame(REAL8TimeSeries *timeData, CHAR *frameType, CHAR *channel, CHAR *ifo)
{
    CHAR fname[100];
    INT4 duration;
    INT8 detectorFlags;
    LALFrameH *frame;
    
    int gpsStart = timeData->epoch.gpsSeconds;
    int gpsEnd = gpsStart + (int)timeData->data->length*timeData->deltaT;
    
    
    /* set detector flags */
    if ( strncmp( ifo, "H2", 2 ) == 0 )
        detectorFlags = LAL_LHO_2K_DETECTOR_BIT;
    else if ( strncmp( ifo, "H1", 2 ) == 0 )
        detectorFlags = LAL_LHO_4K_DETECTOR_BIT;
    else if ( strncmp( ifo, "L1", 2 ) == 0 )
        detectorFlags = LAL_LLO_4K_DETECTOR_BIT;
    else if ( strncmp( ifo, "G1", 2 ) == 0 )
        detectorFlags = LAL_GEO_600_DETECTOR_BIT;
    else if ( strncmp( ifo, "V1", 2 ) == 0 )
        detectorFlags = LAL_VIRGO_DETECTOR_BIT;
    else if ( strncmp( ifo, "T1", 2 ) == 0 )
        detectorFlags = LAL_TAMA_300_DETECTOR_BIT;
    else
    {
        fprintf( stderr, "ERROR: Unrecognised IFO: '%s'\n", ifo );
        exit( 1 );
    }
    
    /* get frame filename */
    duration = gpsEnd - gpsStart;
    snprintf( fname, FILENAME_MAX, "%c-%s-%d-%d.gwf", ifo[0], frameType, gpsStart, duration );
    
    /* set the channel name */
    strncpy(timeData->name, channel, LALNameLength);
    
    /* define frame */
    frame = XLALFrameNew( &timeData->epoch, duration, "LIGO", 0, 1,
                         detectorFlags );
    
    /* add channel to frame */
    //XLALFrameAddREAL8TimeSeriesSimData( frame, timeData );
    XLALFrameAddREAL8TimeSeriesProcData(frame, timeData );
    
    fprintf( stdout, "Writing injection to frame: '%s'\n", fname );
    
    /* write frame */
    if (XLALFrameWrite( frame, fname) != 0)
    {
        fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
        exit( 1 );
    }
    
    /* clear frame */
    XLALFrameFree( frame );
    
    return;
}



