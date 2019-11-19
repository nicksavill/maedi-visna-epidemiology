#include <model.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <omp.h>
#include <regex.h>
#include <sys/types.h>
#include <signal.h>

#define _GNU_SOURCE
 
#define ALLTRT -1
#define ALLIND -1
#define NT2 nt2
#define NR2 (Data->Treatment[nt].nind-1)
#define NTS (ntp == ALLTRT? nt1: ntp)
#define NTE (ntp == ALLTRT? NT2: ntp)
#define NRS (nrp == ALLIND? nr1: nrp)
#define NRE (nrp == ALLIND? NR2: nrp)

#define MAXSTORE MAX(4*Sim->covsamples, 2*Sim->samples)
#define MAXTRTSREPS 1000
#define MAXPARS 100
#define NTEMP 29


enum {NONE, LOW, HIGH};
enum {ALL, FITTED_CONSTANT, FITTED_DERIVED, FITTED_AUXILLARY,
      FITTED_DERIVED_AUXILLARY};
enum {MINTEMP, ADAPTIVE_INDIVIDUAL, ADAPTIVE_TREATMENT, ADAPTIVE_EXPERIMENT,
      NONADAPTIVE};

void Used(int *used, char *line);
void Ignore(int **used, char *line);
void ReadData(Simulation *Sim, struct data *Data);
void Bayes(Simulation *Sim, struct data *Data);
void ConstructLadder(Simulation*, struct sim_individual*, 
		     double, int, int);
void RunMcMC(Simulation* Sim, struct data* Data, int nt, int nr, int phase);
int  InitialiseIndividualCov(Simulation* Sim, struct data* Data, 
			     struct sim_individual *Rep, int nt, int nr, int phase);
void InitialiseSampleCov(Simulation* Sim, struct data *Data);
void MCMC(Simulation* Sim, struct data* Data, int phase, int adjust, 
	  int ntp, int nrp);
void Metropolis(Simulation *Sim, int nt, int nr,
		struct chain *ChainT, double temp, double DlogP,
		struct parameter_block *PB, int phase, struct data *Data, 
		gsl_rng *stream, double *prop_logL, int *change);
double FulllogLikelihood(Simulation *Sim, struct data *Data, double *theta, 
			 int nt, int nr, int *status);
double sumlogPriors(Simulation *Sim, int nt, double *p, double *logP,
		   struct parameterset *Q, double *hyper);
double logPrior(Simulation *Sim, int nt, double p, struct parameter *P, 
		double *hyper);
void Proposal(gsl_rng *stream, double *a, struct parameter_block *PB);
void AdjustCovariance(Simulation *Sim, struct parameter_block *PB, 
		      struct parameterset *Q, int phase); 
void CovMatrix(Simulation *Sim, struct parameter_block *PB);
void ClearThread(struct thread *Thread, Simulation *Sim, struct data *Data, 
		 int ntp, int nrp);
void SetupThread(struct thread *Thread, Simulation *Sim, struct data *Data, 
		 int phase, int ntp, int nrp, int thread);
void SaturatedLikelihood(Simulation *Sim, struct data *Data, int nvar);
int Time_Index(struct data_individual *Ind, int nsteps, double **y, 
	       int *t_index);
int nsteps(Simulation *Sim, struct data_individual *Ind);
void Best(Simulation *Sim, struct data *Data, struct parameterset *Q, 
	  struct chain *ChainT, struct chain *Chain,
	  int ntp, int nrp, int c, int phase);
int Done(Simulation *Sim, struct data *Data, int ntp, int nrp, int phase);
void Default(Simulation *Sim);
double GelmanRubinStats(Simulation *Sim, struct data *Data, FILE *fp);
void ReadParameters(Simulation *Sim, struct parameterset *Q, 
		    int ntrt, int nind, double *p);
void Clear(double **y, int n, int m);
int Use(int type, int use);
void Gibbs(Simulation *Sim, struct data *Data, struct thread *T,
	   gsl_rng *stream);
void SetNAN(double *x, int n);
void PostAnalysis(struct simulation *Sim, struct data *Data);
void Marginals(struct simulation *Sim, struct data *Data);
void PosteriorPrediction(struct simulation *Sim, struct data *Data);
void RecalculateDerivedParameters(Simulation *Sim, struct data *Data);
double Hypervalue(Hyperblock *HB, double *hyper, int pindex, int nt);
void HyperDist(Simulation *Sim, int k, int g, double *hyperlocation, 
	       double *hyperscale);
void SimulatedData(Simulation *Sim, struct data *Data);
void ConstructPiecewiseLine(struct line *L);
void ReturnPiecewiseLine(double x, double *f, double *df, struct line *L);
double NewtonRapherson(double t0, double logr, struct line *F, struct line *G);
void Store(Simulation *Sim, struct chain *Chain, struct chain *ChainT, 
	   struct sim_individual *Ind, int c, struct parameterset *Q, int phase);
int In_set(int x, int *a, int n);
int Check(Simulation *Sim, int k, double logP, double theta, int nt, int nr, 
	  struct parameter *P, double scale);

void SetupHyperBlock(Hyperparameters *HP, struct parameterset *Q, int k, int nt, 
		     int block);

int finite(double);
int omp_get_num_procs(void);
int omp_get_thread_num(void);
long int lrint(double x);
char* strcpy(char*, const char*);
ssize_t getline(char**, size_t*, FILE*);

extern int       choldc_block(double**, long long);
extern int*      integervector(int);
extern int**     integermatrix(int, int);
extern uint**    uintegermatrix(int, int);
extern int***    integermatrix3(int, int, int);
extern void      Int_HeapSort(int*, int);
extern void      HeapSort(double*, int);
extern char*     charvector(int);
extern double    mean(double*, int, int);
extern double*   doublevector(int);
extern double complex*   doublecomplexvector( int );
extern double**  doublematrix(int, int);
extern double*** doublematrix3(int, int, int);
extern double****doublematrix4( int, int, int, int );
extern void      Error(char*);
extern void      ReadParam(char*, void*, char*);
extern void      circmeanvar( double*, int, double*, double* );
extern void      meanvar(double*, int, double*, double*);
extern void      meanvarl( long double*, int, long double*, long double* );
extern void      wmeanvar(long double*, int, long double*, long double*);
/*
extern Matrix *Matrix_Init(Matrix*, char*, int, int);
extern Matrix *Matrix_Destroy(Matrix*);
extern Matrix *Matrix_Copy(Matrix*, Matrix*);
extern Matrix *Matrix_Sub_Copy(Matrix*, Matrix*, int);
extern Matrix *Matrix_Scalar_Multiply(Matrix*, double, Matrix*);
extern Matrix *Matrix_Multiply(Matrix*, Matrix*, Matrix*);
extern Matrix *Matrix_Multiply(Matrix*, Matrix*, Matrix*);
extern Matrix *Matrix_Cholesky(Matrix*, Matrix*);
extern Matrix *Matrix_Invert(Matrix*, Matrix*);
extern Matrix *Matrix_Add(Matrix*, Matrix*, Matrix*);
extern Matrix *Matrix_Add_Diagonal(Matrix*, Matrix*);
extern Matrix *Matrix_Subtract(Matrix*, Matrix*, Matrix*);
extern Matrix *Matrix_Transpose(Matrix*, Matrix*);
extern double Matrix_Log_Determinant(Matrix*, int);
extern void   Matrix_Print(Matrix*);
*/

