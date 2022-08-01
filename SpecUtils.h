#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>

#define sthresh 10.0
#define warm 6.0

#define MM 200000    // iterations of MCMC
#define verbose 0   // set to 0 for quiet run, 1 for verbose
#define dfmin 4.0  // minimum spline spacing
#define dfmax 32.0  // maximum spline spacing
#define smooth 8.0   // moving average
#define tol 0.2     // tolerance for difference in averages
#define linemul 8.0 // how much above the Gaussian floor a line needs to be
#define itype 1   // 0 for amkima splines, 1 for smoothed linear
#define NTHREAD 16  // number of threads
#define ompflag 0 // 1 to use open mp for the likelihood, 0 to not use it
#define Qprint 8.0    // Q used to make output scans


static const gsl_rng_type *rngtype;
static const gsl_rng *rng;


typedef struct {
    int ncof,ioff,joff;
    double *cc,*cr;
} wavefilt;


struct Smooth
{
    int Nsten;
    int Nlook;
    double asmooth;
    double *lookup;
    double dx;
};

void smootset(struct Smooth *smoothline);
void getrangeakima(int k, int Nknot, double *ffit, double Tobs, int *imin, int *imax);
void getrangesmooth(struct Smooth *smoothline, int k, int Nknot, int Ns, double *ffit, double Tobs, int *imin, int *imax);
void setupsline(struct Smooth *smoothline, int Nknot, double *sline, double *ffit, double *slopes, double *bb, double *cc, double *dd);
void delta_setupsline(struct Smooth *smoothline, int iu, int Nknot, double *sline, double *ffit, double *slopes,  double *bb, double *cc, double *dd);
void sline(struct Smooth *smoothline, int istart, int iend, int Nknot, double *farray,  double *SM, double *ffit, double *slinep, double *bb, double *cc, double *dd);
double ltc(struct Smooth *smoothline, double x);
double line(double f, double linef, double lineh, double linew, double deltafmax, double lineQ);

void Inverse(double **M, double **IM, int d);
void whiten(double *data, double *Sn, int N);
void qscan(double *data, double *Sn, double Tobs, int N);
void qscanf(double *data, double *Sn, double Tobs, int N);
void qscanshiftf(double **data, double **Sn, double Tobs, int N);
void qscanmaxf(double *data, double *Sn, double Tobs, int N);
void qscanres(double *data, double *signal, double *Sn, double Tobs, int N);
void tukey(double *data, double alpha, int N);
void tukey_scale(double *s1, double *s2, double alpha, int N);
void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int imin, int imax, int N);
double fourier_nwip(double *a, double *b, double *Sn, int imin, int imax, int N);
void max_array_element(double *max, int *index, double *array, int n);
void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n);
double det(double *A, int N);
double f_nwip(double *a, double *b, int n);
void Ring_all(double *paramx, double *paramy, int d1, int d2, double GMST, gsl_rng * r);
void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m);
void layerC(double *a, double f, double *tf, double *tfR, double *tfI, double Q, double Tobs, double fix, int n);
void SineGaussianC(double *hs, double *sigpar, double Tobs, int N);
void specest(double *data, double *Hf, int N, int Ns, double dt, double fmx, double *SN, double *SM, double *PS);
void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2);
void clean(double *D, double *Draw, double *Hf, double *sqf, double *freqs, double *Sn, double *specD, double *sspecD, double df, double Q, double Tobs, double scale, double alpha, int Nf, int N, int imin, int imax, double *SNR, int pflag);
void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, int N);
void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase);
void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX);
double Getscale(double *freqs, double Q, double Tobs, double fmx, int n, int m);

void makespec(int Nspline, int Nlines, double *ffit, double *Sspline, double *linef, double *lineh,  double *lineQ, double *linew, double *deltafmax, double *SN, double Tobs,  int N);


int *int_vector(int N);
void free_int_vector(int *v);
double *double_vector(int N);
void free_double_vector(double *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);
int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);
double ****double_quad(int N, int M, int L, int K);
void free_double_quad(double ****t, int N, int M, int L);

