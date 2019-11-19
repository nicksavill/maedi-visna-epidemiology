#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <stdbool.h>
#include <complex.h>

/* don't check ranges in Matrix operations and make Matrix call inline */
#define GSL_RANGE_CHECK_OFF
//#define HAVE_INLINE

#define FALSE 0
#define TRUE 1

#ifndef M_PI
#define M_PI            3.141592653589793238462643383280        /* pi */
#endif

#define LN2PI 1.837877066409345339

#define ODEdriver gsl_odeiv2_driver_alloc_y_new
#define rkf45 gsl_odeiv2_step_rkf45
#define rk8pd gsl_odeiv2_step_rk8pd
#define rkck  gsl_odeiv2_step_rkck
#define rk2   gsl_odeiv2_step_rk2

#define EOL -1
       
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define V(x) (x)*(x)

#define FIT_TIME_POINT(x) (x->mode == FIT)
#define RECORD_TIME_POINT(x, u) (x->mode == OUTPUT || (x->mode == FIT && u == x->t[x->i]))
#define NEXT_TIME_POINT(x) if (x->mode == FIT) x->i++
#define RESET_TIME_POINT(x) x->i = 0

enum {FIT, OUTPUT};
  
enum {SUCCESS=0, FAIL=1, ERROR1=1, ERROR2, ERROR3, ERROR4, ERROR5, ERROR6,
      ERROR7, ERROR8, ERROR9, ERROR10};

typedef uint boolean;

// for internal use
// type of parameter
enum {NOTUSED=0, DERIVED=21, FITTED=22, CONSTANT=23, AUXILLARY=24, INDIVIDUAL_CONSTANT=28};
// whether a fitted or auxillary parameter has a hyperparameter associated with it
enum {HYPER=26};
// prior of a fitted parameter that does not have a hyperparameters
enum {NOTAPPLICABLE=0, NORMAL=1, GAMMA=2, SICHI2=3, UNIFORM=4, PARETO=5,
      TRUNCNORMAL=6, EXPONENTIAL=7, BINOMIAL=8, POISSON=9, INVGAMMA=10, 
      UNIFORM0=11, BETA=12, NORMALMEAN=13, CATEGORICAL=14, GEOMETRIC=15,
      USERDEFINED=16, AUXILLARY_HYPER=27, CIRCULAR=29};
// number type of the parameter
enum {REAL=17, INTEGER=18, CATEGORY=19, INT = INTEGER};
// type of data
enum {NOT_AVAILABLE=0, ABOVE_DETECTION=1, BELOW_DETECTION=2, MEASURED=3};

// for use in model files
// prior or type of the declared parameter
enum {Normal        =NORMAL, 
      Gamma         =GAMMA, 
      Sichi2        =SICHI2, 
      Uniform       =UNIFORM,
      Circular      =CIRCULAR,
      Pareto        =PARETO, 
      TruncNormal   =TRUNCNORMAL, 
      Exponential   =EXPONENTIAL, 
      Binomial      =BINOMIAL, 
      Poisson       =POISSON, 
      InvGamma      =INVGAMMA, 
      Uniform0      =UNIFORM0, 
      Beta          =BETA, 
      NormalMean    =NORMALMEAN,
      Categorical   =CATEGORICAL, 
      Geometric     =GEOMETRIC, 
      UserDefined   =USERDEFINED, 
      HyperP         =HYPER, 
      Constant      =CONSTANT, 
      Auxillary     =AUXILLARY, 
      Auxillary_Hyper =AUXILLARY_HYPER, 
      Derived       =DERIVED,
      Individual_Constant=INDIVIDUAL_CONSTANT};    
enum {Nor           =NORMAL, 	  
      Gam 	    =GAMMA, 	  
      Sic 	    =SICHI2, 	  
      Unif 	    =UNIFORM, 	  
      Circ          =CIRCULAR,
      Pare 	    =PARETO, 	  
      TNorm 	    =TRUNCNORMAL, 
      Exp 	    =EXPONENTIAL, 
      Bin 	    =BINOMIAL, 	  
      Poi 	    =POISSON, 	  
      IGam 	    =INVGAMMA, 	  
      Unif0 	    =UNIFORM0, 	  
      Be 	    =BETA, 	  
      NorM 	    =NORMALMEAN,  
      Cat 	    =CATEGORICAL, 
      Geom 	    =GEOMETRIC,   
      UD 	    =USERDEFINED, 
      Hyp 	    =HYPER, 	  
      Cst 	    =CONSTANT, 	  
      Aux 	    =AUXILLARY,   
      Aux_Hyp       =AUXILLARY_HYPER, 
      Der	    =DERIVED,
      IndCst        =INDIVIDUAL_CONSTANT};    

// number type of the parameter
enum {Real=REAL, Integer=INTEGER, Category=CATEGORY, Int=INTEGER};
// type of data
enum {NA=0, Above_Detection=1, Below_Detection=2, Measured=3};

typedef struct {
  uint mode, n, i;
  double **t;
} TimePoints;

typedef int (*function_ptr)(int, int, double*, double*, double*, TimePoints*);

typedef struct {
  int n_func;
  function_ptr *function_list;
} functions_t;

typedef struct {
  uint n_timepoints;
  double *time;
  double **value;
} Data_t;

struct line {
  int n;
  double *x;
  double *y;
  double *m;
  double *c;
};

struct data_individual {
  int n;
  int np;
  uint **value;
  double **Y;
  double saturated;
};

struct data_treatment {
  int nind;
  int *usevar;
  char levels[1000];
  struct data_individual *Individual;
};

struct data {
  int ntrt;
  int nvar;
  int maxnind;
  int maxnp;
  char filename[100];
  struct data_treatment *Treatment;
};

typedef struct treatmentlist {
  int n;
  int *treatment;
} Treatments;

typedef struct hyperblock {
  /* set in CreateBlock */
  int use;

  /* set in AddHyperToTreatment */
  int nhypers_used;
  int *hyper_used;
  int **nhypers_with_treatment;
  int ***hyper_with_treatment;
  int nassociated_pars;
  int *associated_par;
  int *ntreatments_with_parameter;
  int **treatment_with_parameter;

  /* set in Bayes */
  int hdistribution;
  int pdistribution;
  int *dependent_block;
} Hyperblock;

struct hyper {
  int index;
  int use;
  int parno;
  int hdist;
  int number;
  char name[100];
  double location;
  double scale;
  double shape;
};

typedef struct allhypers {
  int init;
  int ntrt;
  int *nind;
  int *trtused;
  int nhyper_parameters;
  struct hyper *H;
  int nblocks;
  Hyperblock *block;
  double *hyper;         //used for initialisation of thread hyper values
  double **hyper_store;
} Hyperparameters;

struct variable {
  int index;
  char name[100];
};
  
struct variables {
  int    variables;
  double scaletime;
  struct variable *W;
};

struct parameter {
  int index;
  int number;
  int distribution;
  int circular;
  int type;
  int hashyper;
  int hyper_location_block_index;
  int hyper_shape_block_index;
  int hyper_scale_block_index;
  int ncategories;
  double *p; // used for categorical parameters
  double location;
  double scale;
  double shape;
  char name[100];
  double (*(*link))(double);
};
  
struct blocks {
  int nblocks;
  int **pblock;
};

struct parameterset {
  int    pass;
  int    parameters;   // maximum parameter index
  int    npar;         // number of fitted parameters
  int    ninteger;
  int    *integer;
  int    nauxillary;
  int    *auxillary;
  struct parameter *P;
  struct blocks Blocks;
};

struct parameter_block {
  int categorical;
  int done;
  int npar;
  int *pblock;
  int covfinal;
  int covstart;
  int inrange;
  int storeno;
  int error;
  double accept;      /*mcmc acceptances*/
  double proposed;    /*mcmc iterations*/
  double covscale;    /*covariance matrix scaling factor*/
  //  double **L;         /*covariance matrix*/
  gsl_matrix *L;
  double **theta_store;
};

struct chain {
  int done;       /*if current individual is ready for next phase*/
  int loop;
  /* parameters for threads */
  double logmaxL;
  double logL;
  double *logP;
  double logpost;
  double *theta;       /*current parameter values*/
  double minlogL;
  double prop_logL;
  double prop_logP;
  double prop_logpost;
  int storeno;
  int laststoreno;
  double accept;      /*mcmc acceptances*/
  double proposed;    /*mcmc iterations*/
  double **theta_store;
  double *logLstore;
  double *logpoststore;
  double exchangeaccept;      /*mcmc acceptances*/
  double exchangeproposed;    /*mcmc iterations*/
  struct parameter_block *PB;
};

struct sim_individual {
  int done;
  int chains;
  int *chainno;
  double *temp;
  struct chain *Chain;
  struct parameterset Q;
};

typedef struct {
  int n;
  int m;
  int **a;
} Matrix;

struct sim_treatment {
  struct sim_individual *Individual;
};

struct thread {
  double *hyper;
  struct sim_treatment *Treatment;
  gsl_rng *stream;
};

typedef struct simulation {
  Hyperparameters HP;
  int *trtused;
  int **indused;
  int ntp;            // individual with most parameters
  int nrp;            // individual with most parameters
  int rdp;
  int sim;
  int maxc;
  int seed;
  int thin;
  int mcmc;
  int popd;
  int marg;
  int maxnp;
  int error;
  int range;
  int chains;
  int notime;
  int copyind;
  int threads;
  int verbose; /* 0: no output, 1: output chain 0+lowest exchangerate, 
		  2: all chains*/
  int derived;
  int discard;
  int samples;
  int storeno;
  int gsl_error;
  int population;
  int covsamples;
  int growthcurves;
  int parsperblock;
  int phase1adjust;
  int phase2adjust;
  int mincovsamples;
  int ignorewarning;
  int usehypermeans;
  int sample_parameters;
  char expno[50];
  char outexpno[50];
  char *parsfile;
  char *tempfile;
  char logfile[1000];
  char factors[1000];
  double xrate;
  double timeend;
  double lowrate;
  double highrate;
  double minprior;
  double initscale;
  double covthresh;
  double timestart;
  double scalefactor;
  double saturatedlogL;
  double minexchangerate;
  struct sim_individual Sample;
  struct sim_treatment *Treatment;
  struct variables V;
} Simulation;

extern char expno[];

extern bool **boolmatrix( int, int );
extern int* integervector( int );
extern uint* uintegervector( int );
extern uint** uintegermatrix( int, int );
extern double* doublevector( int );
extern double complex* doublecomplexvector( int );
extern bool* boolvector( int );
extern double** doublematrix( int, int );
extern double*** doublematrix3( int, int, int );
extern double**** doublematrix4( int, int, int, int );
extern char* charvector( int );
extern char** charmatrix( int, int );
extern double factln( int );
extern double gammln( double );
extern void binorm( double*, double, double, double, double, double );
extern long int lrint(double);
extern double pow10(double);
extern double round(double);
extern double rtnewt(void(*)(double, double*, double*, double*), 
		     double, double, double, double*);
extern double rtbis( double (*)(double, double*), double, double, 
		     double, int*, double*);

extern void *malloc(size_t);
extern void exit(int);
extern void free(void*);
extern void Error(char*);
extern void ReadParam( char*, void*, char* ) ;

extern uint   binomialRNG(uint, double, const gsl_rng*);
extern void   multinomialRNG(uint, const double*, size_t, unsigned int*, const gsl_rng*);
extern int    poissonRNG(double, const gsl_rng*);
extern int    uniform_intRNG(int, int, const gsl_rng*);
extern double normalRNG(double, double, const gsl_rng*);
extern double betaRNG(double, double, const gsl_rng*);
extern double chisqRNG(double, const gsl_rng*);
extern double scaleinvchisqRNG(double, double, const gsl_rng*);
extern double invgammaRNG(double, double, const gsl_rng*);
extern double gammaRNG(double, double, const gsl_rng*);
extern double uniformRNG(double, double, const gsl_rng*);
extern double paretoRNG(double, double, const gsl_rng*);
extern void   dirichlet_multinomialRNG(uint, const double*, size_t, uint*, const gsl_rng*);
extern void   dirichletRNG(const double*, size_t, double*, const gsl_rng*);
extern uint   beta_binomialRNG(uint, double, double, const gsl_rng*);
extern uint   negative_binomialRNG(double, double, const gsl_rng*);
extern uint   beta_negative_binomialRNG(uint, double, double, const gsl_rng*);
/* Matrix *wishartRNG(int, int, Matrix*, Matrix*, const gsl_rng*); */
/* void   MVNRNG(int, Matrix*, Matrix*, Matrix*, const gsl_rng*); */

extern double poissonCMF(const double, const double, int);
extern double gamma_scalePDF(const double, const double, const double);
extern double betaPDF(const double, const double, const double);
extern double exponentialPDF(const double, const double);
extern double poissonPMF(const double, const double);
extern double geometricPDF(const double, const double);
extern double trnormalPDF(const double, const double, const double);
extern double normalCDF(const double, const double, const double, int);
extern double normalPDF(const double, const double, const double);
extern double binomialPMF(const double, const double, const double);
extern double multinomialPMF(const unsigned int, const double*, const double*);
extern double uniformPDF(const double, const double, const double);
extern double logisticPDF(const double, const double, const double);
extern double betabinomialPMF(const double, const double, const double, const double);
//extern double MVNPDF(Matrix*, Matrix*, Matrix*, const double, int);
extern double von_mises_cdf ( double x, double a, double b );
extern double von_mises_cdf_inv ( double cdf, double a, double b );

extern Matrix *MatrixInit(Matrix*, int, int, int*);
extern Matrix *MatrixDestroy(Matrix*);
extern double identity(double);
extern double logit(double);
extern double invlogit(double);


void function_list(functions_t*);
// int    function(int, int, double*, double*, double*, TimePoints*);
int    Model(const uint, const uint, const uint, double*, double**,
	     TimePoints*);
//int    Model(int, int, int, double*, double**);
void   Block(struct parameterset*, int, ...);
void   Par(struct parameterset*, int, int, char*, int, ...);
void   Create(Hyperparameters *H, int *h, int block, int par,
	     Matrix M, char *name, int type, double p1, double p2);
void   Hyper(Hyperparameters *HP, int index, int number, char *name, int distribution,
	     double shape_or_loc, double scale);
void   CreateHyperBlock(Hyperparameters *HP, int index);
void   AddHyperToATreatment(Hyperparameters *HP, int index, int h, int pindex,
			    int treatment);
void   AddHyperToTreatments(Hyperparameters *HP, int index, int h, int pindex, ...);
void   HyperParameters(Treatments, Hyperparameters*);
void   Var(struct variables*, int, char*); 
void   GlobalParameters(void);
void   Variables(struct variables*);
void   OutputModel(int, int, double*, double*);
void   OutputData(int, int, double*, double*, double*, uint*, int);
void   DerivedHyperParameters(int trt, double *p);
void   DerivedParameters(int, int, int, double*, double**, TimePoints*);
void   Parameters(int, int, double*, struct parameterset*); 
void   PredictData(int, int, double*, double*, double*, gsl_rng*);
void   SimulateData(int, int, double*, double*, double*, double*, uint*, gsl_rng*);
void SaturatedModel(int trt, int ind, double **V, double *p, 
        const TimePoints *Data);
// void   SaturatedModel(int, int, double*, double*, double*, uint*); 
void   Residual(int, int, double*, double*, double*, double*, uint*);
void   Equations(double*, double*, int, double*, double*, int, double);
double timestep(void);
double logLikelihood(const uint, const uint, const double*,
		     double**, const TimePoints*);
// double logLikelihoodI(const uint trt, const uint ind, const double *p, 
// 		      double **V, double *Data, const uint *value,
		      // const uint N);
void WAIC(int, int, double**, double**, double*, const TimePoints*);
// void WAIC(int trt, int ind, double *lnL, double *V,
// 	  double *p, double *Data, uint *value);
//double logLikelihood(int, int, double*, double*, double*, int*);
double UserPDF(double);
void   ScaleTime(struct variables*, double);
