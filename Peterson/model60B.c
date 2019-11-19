#include <model.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#define P_inv_gamma gsl_cdf_gamma_Pinv // invverse cumulative gamma dist
#define P_inc_gamma gsl_sf_gamma_inc_P // regularised incomplete Gamma function
#define lngamma gsl_sf_lngamma       // ln Gamma function

#define V0 9
#define NUMBER_SHEEP 28
#define NUM_DONORS 9
#define SHEEP_FIRST_TEST_WEEK 1

// month defined as 365/12 days 
// week defined as 7 days
# define WEEKS_PER_MONTH ((365/12)/7) 

#include "Blessumer_clean_data.h"

typedef struct {
  double a, b;
  double *lambda;
  double *Lambda;
  double **W4;
  double *W2;
} parameters_t;

double u(int, double*);
double lngexpan(double, double, int*);
double ln_inc_gamma_star(double, double);
double w(int, int, int, parameters_t*);
double Prob_seroconvert(int, parameters_t*);
extern double R0(double, double, double, int, double*);

void Variables(struct variables *V) 
{
  Var(V, 0, "week");
  Var(V, 1, "field infected");
  Var(V, 2, "new positives");
  Var(V, 3, "cumulative new positives");
  Var(V, 4, "seroconversion rate");
  Var(V, 5, "number infectious");
  Var(V, 6, "number infectious");
  Var(V, 6, "housed infected");
  Var(V, 7, "infected");
  Var(V, 8, "infected");
  Var(V, V0, "log_p_pos"); 
}

void Parameters(int trt, int ind, double *p, struct parameterset *Q) 
{
  Par(Q,  0, Real, "beta_{field}", Exponential, 0.0008);
  Par(Q,  1, Real, "beta_{housed}", Gamma, 2.0, 0.1);
  Par(Q,  2, Real, "a", Exponential, 3.2);
  Par(Q,  3, Real, "k1", Gamma, 2.0, 4.0);
  // Par(Q,  0, Real, "beta_{field}", Exponential, 0.0082/WEEKS_PER_MONTH);
  // Par(Q,  1, Real, "beta_{housed}", Gamma, 2.0, 0.01);
  // Par(Q,  2, Real, "mu", Gamma, 2.0, 1.5);
  // Par(Q,  3, Real, "sigma", Gamma, 2.0, 1.5);

  // Par(Q,  8, Real, "a", Der);
  // Par(Q,  9, Real, "b", Der);
  // Par(Q, 10, Real, "beta_{housed}/beta_{field}", Der);
  // Par(Q, 12, Real, "E[inf_{field}]", Der);
  // Par(Q, 13, Real, "E[inf_{housed}]", Der);
  // Par(Q, 14, Real, "E[inf]", Der);
  // Par(Q, 15, Real, "s_{low}", Der);
  // Par(Q, 16, Real, "s_{high}", Der);
}

int Model(const uint trt, const uint ind, const uint nweeks, double *p, double **V, TimePoints *Data)
{
  int i, j, m, *num_infected;
  double **p_seropos = V, *Einf_field, *Einf_housed;
  parameters_t pars;

  double beta_field = p[0]/WEEKS_PER_MONTH;
  double beta_housed = p[1]/WEEKS_PER_MONTH;
  double S_a = p[2];
  double S_b = S_a/p[3]/WEEKS_PER_MONTH;
  int seroconversion = lrint(p[3]);
  
  for (m = 0; m < nweeks; m++)
    V[m][0] = m;

  int *num_infectious = integervector(nweeks);
  pars.Lambda = doublevector(nweeks);
  pars.lambda = doublevector(nweeks);
  pars.W4 = doublematrix(nweeks, nweeks);
  pars.W2 = doublevector(nweeks);
  pars.a = S_a;
  pars.b = S_b;

  // ewe to ewe transmission
  // the 10 donor sheep are infectious at start of experiment
  for (j = 0; j < NUM_DONORS; j++)
    for (m = 0; m <= sheep_last_week[j]; m++)
        num_infectious[m]++;

  for (m = 0; m < nweeks; m++)
    if (m == 18 || m == 19)
      pars.lambda[m] = num_infectious[m]*beta_housed;
    else
      pars.lambda[m] = num_infectious[m]*beta_field;

  // cumulative sum of lambda
  pars.Lambda[0] = 0;
  for (m = 1; m < nweeks; m++)
    pars.Lambda[m] = pars.Lambda[m-1]+pars.lambda[m-1];

  for (m = 0; m < nweeks; m++) {
    for (i = 0; i < m; i++)
      pars.W4[i][m] = w(i, i+1, m+1, &pars) 
                    - w(i, i+1, m,   &pars) 
                    - w(i, i,   m+1, &pars)  
                    + w(i, i,   m,   &pars);
    pars.W2[m] = w(m, m, m+1, &pars);
  }
 
  // for each testing date calculate the probability that a sheep tests positive since the last test
  for (m = 0; m < nweeks; m++)
    p_seropos[m][V0] = Prob_seroconvert(m, &pars);



  /******************** OUTPUT ***********************/
  if (Data->mode == OUTPUT) {
    for (m = 0; m < nweeks; m++)
       V[m][0] = m;

    // expected number of new positives on each test date
    for (j = NUM_DONORS; j < NUMBER_SHEEP; j++)
      for (i = 0; i < sheep_last_test_week[j]-1; i++)
        V[i+1][2] += p_seropos[i][V0];

    // cumulative expected number of new positives
    V[0][3] = V[0][2];
    for (i = 0; i < nweeks-1; i++)
      V[i+1][3] = V[i][3]+V[i+1][2];

    // seroconversion rate (new positives per week)
    // V[0][4] = 0;
    // for (i = 1; i < Data->n; i++)
    //   V[i][4] = V[i][2]/(test_weeks[i]-test_weeks[i-1]);

    // // number of infectious animals
    // for (i = 0; i < Data->n; i++)
    //   V[i][5] = num_infectious[test_weeks[i]];

    // Einf_field = doublevector(nweeks);
    // Einf_housed = doublevector(nweeks);
    // for (j = NUM_DONORS; j < NUMBER_SHEEP; j++)
    //   for (m = 0; m <= SHEEP_LAST_WEEK; m++)
    //     if (test_weeks[m] == 19 || test_weeks[m] == 20)
    //       Einf_housed[m] += (1.-exp(-pars.lambda[m]))*u(m, pars.Lambda);
    //     else
    //       Einf_field[m] += (1.-exp(-pars.lambda[m]))*u(m, pars.Lambda);

    // for (m = 1; m < nweeks; m++) {
    //   Einf_field[m] += Einf_field[m-1];
    //   Einf_housed[m] += Einf_housed[m-1];
    // }
    // for (i = 0; i < Data->n; i++) {
    //   V[i][1] = Einf_field[test_weeks[i]];
    //   V[i][6] = Einf_housed[test_weeks[i]];
    //   V[i][7] = V[i][1]+V[i][6];
    // }
    // free(Einf_field);
    // free(Einf_housed);

    // number infected
    // num_infected = integervector(nweeks);
    // for (j = 0; j < NUMBER_SHEEP; j++)
    //   if (sheep_pos_week[j] != -1)
    //     for (m = max(0, sheep_pos_week[j]-seroconversion); m <= SHEEP_LAST_WEEK; m++)
    //       num_infected[m]++;

    // for (i = 0; i < Data->n; i++)
    //   V[i][8] = num_infected[test_weeks[i]];
    // free(num_infected);
  }

  free(num_infectious);
  free(pars.lambda);
  free(pars.Lambda);
  free(pars.W4[0]); free(pars.W4);
  free(pars.W2);
  return SUCCESS;
}

double P_pos(int s1, int s2, double **V)
{
  int m;
  double sum = 0;
  for (m = s1; m < s2; m++)
    sum += V[m][V0];
  return sum;
}

double logLikelihood(const uint trt, const uint ind, const double *p, double **V, const TimePoints *Data)
{
  int j, pos_idx;
  double lnL = 0;

  for (j = NUM_DONORS; j < NUMBER_SHEEP; j++)
    if (sheep_pos_index[j] != 0) {
      pos_idx = sheep_pos_index[j];
      if (pos_idx > 0)
        // ewe is first seropositive at start of month m, therefore it seroconverted
        // sometime in the months from the last test to the month prior to seroconverting
        lnL += log(P_pos(lrint(Data->t[pos_idx-1][0]), sheep_pos_week[j], V));
      else
        // ewe never seroconverted
        lnL += log(1. - P_pos(SHEEP_FIRST_TEST_WEEK, sheep_last_test_week[j], V));
    }
  return lnL;
}

void OutputModel(int trt, int ind, double *Output, double *V)
{
  int i;
  for (i = 1; i < V0; i++)
    Output[i] = V[i];
}

void OutputData(int trt, int ind, double *Output, double *Var, double *Data, uint *value, int index)
{
  Output[2] = Data[29];
  Output[3] = Data[30];
  Output[4] = Data[31];
  Output[5] = Data[32];
}

int f0(int trt, int ind, double *x, double *y, double *p, TimePoints *TP)
{
  // R0 2 month housed
  int i, n = 2;
  if (x == NULL) return n;
  for (i = 0; i < n; i++) {
    x[i] = 100*i+1;
    y[i] = R0(x[i], 2., 0.1, 5, p);
  }
  return n;
}

int f1(int trt, int ind, double *x, double *y, double *p, TimePoints *TP)
{
  // R0 0 months housed
  int i, n = 2;
  if (x == NULL) return n;
  for (i = 0; i < n; i++) {
    x[i] = 100*i+1;
    y[i] = R0(x[i], 0., 0.1, 5, p);
  }
  return n;
}

int f2(int trt, int ind, double *x, double *y, double *p, TimePoints *TP)
{
  // R0 for 25 sheep housed for different durations
  int i, n = 14;
  if (x == NULL) return n;
  for (i = 0; i < n; i++) {
    x[i] = i;
    y[i] = R0(25., x[i]*12./365., 0.1, 5, p);
  }
  return n;
}

int f3(int trt, int ind, double *x, double *y, double *p, TimePoints *TP) 
{
  // CDF of normal seroconversion period
  int i, n = 1000;
  if (x == NULL) return n;
  double S_a = p[2];
  double S_k1 = p[3];

  for (i = 0; i < n; i++) {
    x[i] = (double)i/20.;
    // y[i] = gsl_ran_gamma_pdf(x[i], a, 1./b);
    y[i] = gsl_cdf_gamma_P(x[i], S_a, S_k1);
  }
  return n;
}

void function_list(functions_t *F)
{
  F->n_func = 0;
  // F->function_list = (function_ptr*) malloc(F->n_func*sizeof(function_ptr));
  // F->function_list[0] = f0;
  // F->function_list[1] = f1;
  // F->function_list[2] = f2;
  // F->function_list[3] = f3;
}

double u(int m, double *Lambda)
{
  // The probability a sheep is uninfected at the start of week m
  return exp(-Lambda[m]);
}

double lngseries(double a, double x, int *ierr)
{
  double giant = HUGE_VAL/1000., eps = 1e-15;
  double t = 1./a, v = t;
  double p, q, lnigam;
  int k = 0;
  *ierr = 0;

  while ((fabs(t/v) > eps) && *ierr == 0) {
    p = (a+k)/(a+k+1);
    q = k+1;
    t *= -x*p/q;
    v += t;
    k += 1;
    if (t > giant)
      *ierr = 1;
  }
  if (*ierr == 0)
    if (lngamma(a) < log(giant))
      lnigam = log(v)-lngamma(a);
    else {
      lnigam = 0;
      *ierr = 1;
    }
  else {
    lnigam = 0;
    *ierr = 1;
  }
  return lnigam;
}

double lngexpan(double a, double x, int *ierr)
{
  double giant = HUGE_VAL/1000., eps = 1e-15;
  double t = 1, v = t;
  double p, lnigam;
  int k = 0;
  *ierr = 0;

  if (x > log(giant)) {
      *ierr = 1;
      return 0;
  }
  while ((fabs(t/v) > eps) && *ierr == 0) {
    p = 1-a+k;
    t *= p/x;
    v += t;
    k += 1;
    if (t > giant)
      *ierr = 1;
  }
  if (*ierr == 0)
    if (lngamma(a) < log(giant))
      lnigam = log(v)+x-log(x)-lngamma(a);
    else {
      lnigam = 0;
      *ierr = 1;
    }
  else {
    lnigam = 0;
    *ierr = 1;
  }
  return lnigam;
}

double ln_inc_gamma_star(double a, double z)
{
  int ierr;
  double lnigam;

  if (z < -50) {
    lnigam = lngexpan(a, -z, &ierr);
    if (ierr != 0)
      printf("error1\n");
  } 
  else {
    lnigam = lngseries(a, z, &ierr);
    if (ierr != 0)
      printf("error2\n");
  }
  return lnigam;
}

double w(int n, int t, int s, parameters_t *p)
{
  double a = p->a;
  double l = p->lambda[n];
  double b = p->b;
  double r = (b-l)/b;
  double Z, z, v;

  Z = b*(s-t);
  if (Z == 0)
    return 0;

  z = r*Z;
  if (z <= 0)
    v = exp(a*log(Z) - l*Z/b + ln_inc_gamma_star(a, z));
  else
    v = exp(-a*log(r) - l*Z/b + log(P_inc_gamma(a, z)));

  return exp(-l*(t-n))*(v - P_inc_gamma(a, Z));
}

double Prob_seroconvert(int m, parameters_t *p)
{
  int n;
  double P_m = 0;

  for (n = 0; n < m; n++)
    P_m += u(n, p->Lambda)*p->W4[n][m];
  P_m -= u(m, p->Lambda)*p->W2[m];

  return min(max(0, P_m), 1);
}

void DerivedParameters(int trt, int ind, int nmonths, double *p, double **V, TimePoints *Data) 
{
  // int i, j, m;
  // double a, b;
  // double beta_field = p[0];
  // double beta_housed = p[1];
  // double mu = p[2];
  // double sd = p[3];
  // int seroconversion = lrint(mu);
  // int latency = lrint(p[4]);
  // int nweeks = lrint(Data->t[Data->n-1][0])+1;
  // int *test_weeks = integervector(Data->n);
  // int *num_infectious = integervector(nweeks);
  // double *Lambda = doublevector(nweeks);
  // double *lambda = doublevector(nweeks);

  // p[8]  = a = pow(mu/sd, 2.);
  // p[9]  = b = a/mu;
  // p[10] = beta_housed/beta_field;
  // p[15] = P_inv_gamma(0.025, a, 1./b);
  // p[16] = P_inv_gamma(0.975, a, 1./b);

  // for (i = 0; i < Data->n; i++)
  //   test_weeks[i] = lrint(Data->t[i][0]);

  // for (j = 0; j < NUM_DONORS; j++)
  //   for (m = 0; m <= SHEEP_LAST_WEEK; m++)
  //       num_infectious[m]++;

  // for (j = NUM_DONORS; j < NUMBER_SHEEP; j++)
  //   if (sheep_pos_week[j] > 0)
  //     for (m = max(0, sheep_pos_week[j]-seroconversion+latency); m <= SHEEP_LAST_WEEK; m++)
  //       num_infectious[m]++;

  // for (m = 0; m < nweeks; m++)
  //   if (test_weeks[m] == 19 || test_weeks[m] == 20)
  //     lambda[m] = num_infectious[m]*beta_housed;
  //   else
  //     lambda[m] = num_infectious[m]*beta_field;

  // Lambda[0] = 0;
  // for (m = 1; m < nweeks; m++)
  //   Lambda[m] = Lambda[m-1]+lambda[m-1];

  // p[12] = 0;
  // p[13] = 0;
  // for (j = NUM_DONORS; j < NUMBER_SHEEP; j++)
  //   for (m = 0; m <= SHEEP_LAST_WEEK; m++)
  //     if (test_weeks[m] == 19 || test_weeks[m] == 20)
  //       p[13] += (1.-exp(-lambda[m]))*u(m, Lambda);
  //     else
  //       p[12] += (1.-exp(-lambda[m]))*u(m, Lambda);
  // p[14] = p[12]+p[13];
  
  // free(lambda);
  // free(Lambda);
  // free(num_infectious);
  // free(test_weeks);
}

void WAIC(int trt, int ind, double **lnL, double **V, double *p, const TimePoints *Data) {}
void PredictData(int trt, int ind, double *Output, double *V, double *p, gsl_rng *stream) {}
void SimulateData(int trt, int ind, double *Output, double *V, double *p, double *Data, uint *value, gsl_rng *stream) {}
double timestep(void) {return 1;}
void GlobalParameters() {}
double UserPDF(double x) {return 0;}
void HyperParameters(Treatments T, Hyperparameters *H) {}
void SaturatedModel(int trt, int ind, double **V, double *p, const TimePoints *Data) {}
void Residual(int trt, int ind, double *R, double *V, double *p, double *Data, uint *value) {}
