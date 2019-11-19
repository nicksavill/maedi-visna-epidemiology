#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <stdbool.h>

#define FALSE 0
#define TRUE 1

#ifndef M_PI
#define M_PI            3.141592653589793238462643383280        /* pi */
#endif

#define LN2PI 1.83787706640934548356065947

extern char parameterfilename[50];
extern char parameterlogfilename[50];

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
enum {NOT_AVAILABLE=0, ABOVE_DETECTION=1, BELOW_DETECTION=2, MEASURED=3};

typedef struct {
  int n;
  int m;
  int **a;
} Matrix;

void ReadParam( char*, void*, char* ) ;
void ExperimentNumber( char*, char* );
void Error( char* );
void FileError( char*, char* );
double ***doublematrix3( int, int, int );
double **doublematrix( int, int );
bool **boolmatrix( int, int );
int **integermatrix( int, int );
uint **uintegermatrix( int, int );
bool *boolvector( int );
uint *uintegervector( int );
char **charmatrix( int, int );
double *doublevector( int );
double complex *doublecomplexvector( int );
int *integervector( int );
char *charvector( int );

void Int_HeapSort(int*, int);
void HeapSort( double*, int );
void circmeanvar( double*, int, double*, double* );
void meanvar( double*, int, double*, double* );
void meanvarl( long double*, int, long double*, long double* );
void wmeanvar( long double *x, int n, long double *m, long double *s2 );
double mean( double*, int, int );
int choldc_block(double**, long long);

uint   binomialRNG(uint, double, const gsl_rng*);
void   multinomialRNG(uint, const double*, size_t, unsigned int*, 
		      const gsl_rng*);
int    poissonRNG(double, const gsl_rng*);
int    uniform_intRNG(int, int, const gsl_rng*);
double normalRNG(double, double, const gsl_rng*);
double betaRNG(double, double, const gsl_rng*);
double chisqRNG(double, const gsl_rng*);
double scaleinvchisqRNG(double, double, const gsl_rng*);
double invgammaRNG(double, double, const gsl_rng*);
double gammaRNG(double, double, const gsl_rng*);
double uniformRNG(double, double, const gsl_rng*);
double paretoRNG(double, double, const gsl_rng*);
void   dirichlet_multinomialRNG(uint, const double*, size_t, uint*,
				const gsl_rng*);
void   dirichletRNG(const double*, size_t, double*, const gsl_rng*);
uint   beta_binomialRNG(uint, double, double, const gsl_rng*);
uint   negative_binomialRNG(double, double, const gsl_rng*);
uint   beta_negative_binomialRNG(uint, double, double, const gsl_rng*);
/* Matrix *wishartRNG(int, int, Matrix*, Matrix*, const gsl_rng*); */
/* void   MVNRNG(int, Matrix*, Matrix*, Matrix*, const gsl_rng*); */

double poissonCMF(const double, const double, int);
double gamma_scalePDF(const double, const double, const double);
double betaPDF(const double, const double, const double);
double exponentialPDF(const double, const double);
double poissonPMF(const double, const double);
double geometricPDF(const double, const double);
double trnormalPDF(const double, const double, const double);
double normalCDF(const double, const double, const double, int);
double normalPDF(const double, const double, const double);
double binomialPMF(const double, const double, const double);
double multinomialPMF(const unsigned int, const double*, const double*);
double uniformPDF(const double, const double, const double);
double logisticPDF(const double, const double, const double);
double betabinomialPMF(const double, const double, const double, const double);
/* double MVNPDF(Matrix*, Matrix*, Matrix*, const double, int); */

double identity(double);
double logit(double);
double invlogit(double);

Matrix *MatrixInit(Matrix*, int, int, int*);
Matrix *MatrixDestroy(Matrix*);
/* Matrix *Matrix_Init(Matrix*, char*, int, int); */
/* Matrix *Matrix_Fill(Matrix*, int, int, ...); */
/* Matrix *Matrix_Destroy(Matrix*); */
/* Matrix *Matrix_Copy(Matrix*, Matrix*); */
/* Matrix *Matrix_Sub_Copy(Matrix*, Matrix*, int); */
/* Matrix *Matrix_Scalar_Multiply(Matrix*, double, Matrix*); */
/* Matrix *Matrix_Multiply(Matrix*, Matrix*, Matrix*); */
/* Matrix *Matrix_Multiply(Matrix*, Matrix*, Matrix*); */
/* Matrix *Matrix_Cholesky(Matrix*, Matrix*); */
/* Matrix *Matrix_Invert(Matrix*, Matrix*); */
/* Matrix *Matrix_Add(Matrix*, Matrix*, Matrix*); */
/* Matrix *Matrix_Add_Diagonal(Matrix*, Matrix*); */
/* Matrix *Matrix_Subtract(Matrix*, Matrix*, Matrix*); */
/* Matrix *Matrix_Transpose(Matrix*, Matrix*); */
/* double Matrix_Log_Determinant(Matrix*, int); */
/* void   Matrix_Print(Matrix*); */

double rtnewt(void(*)(double, double*, double*, double*), 
	      double, double, double, double*);
double rtbis( double (*)(double, double*), double, double, double, int*,
	      double* );

size_t strlen(const char*);
char *strncpy(char*, const char*, size_t);
char *strcpy(char*, const char*);
char *strcat(char*, const char*);
char *strchr ( const char*, int );
int  strcmp (const char*, const char*);
void *calloc(size_t, size_t);
void *malloc(size_t);
void free(void*);
void exit(int);

double von_mises_cdf ( double x, double a, double b );
double von_mises_cdf_inv ( double cdf, double a, double b );

