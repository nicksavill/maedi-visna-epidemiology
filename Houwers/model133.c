// The first test on 25/3/80 contains no information about
// infection because ewes had not been exposed to the donors at this time
// The first test date is end of March which we consider start of April

// Time is counted in months at the start of each month, so month 0 is 1 April 1980
// A seroconversion observed at the start of month m means the seroconversion 
// occurred before the start of month m, eg in month m-1 or before

// Sheep are assumed to be born on the 1st April of each year, ie months 0, 12 and 24
// A sheep's last month is the last month it was tested. We assume that it is immediately removed
// This means that a +ve sheep is not infectious on its last_month
// We assume, however, that a sheep is potentially infected and infectious at most one month prior to
// its first positive month

// Donor sheep 72 is assumed infectious at the start, because its dam tests positive within 3 months of birth
// Donor 22 is assumed to be latent at the start due to it being removed after 3.5 years,
// so it may have just been infected at the start

#include <model.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "sheep_data.h"

#define P_inv_gamma gsl_cdf_gamma_Pinv // invverse cumulative gamma dist
#define P_inc_gamma gsl_sf_gamma_inc_P // regularised incomplete Gamma function
#define lngamma gsl_sf_lngamma       // ln Gamma function
#define NUMBER_DONORS 2
#define NUMBER_SHEEP 45
#define INTRO_MONTHS 3
#define V0 11

const uint intro_months[] = {0, 12, 24};

typedef struct {
  double a;
  double *lambda;
  double *Lambda;
  double *b;
  double *B;
  double **W4;
  double *W2;
} parameters_t;

double u(int, int, double*);
double lngexpan(double, double, int*);
double ln_inc_gamma_star(double, double);
double w(int, int, int, int, parameters_t*);
double Prob_seroconvert(int, int, parameters_t*);
extern double R0(double, double, double, int, double*);

void Variables(struct variables *V) 
{
  int i;
  Var(V,  0, "month");
  Var(V,  1, "expected number infected");
  Var(V,  2, "expected number infectious");
  Var(V,  3, "cumulative seroconversions");
  Var(V,  4, "seroconversion rate");
  Var(V,  5, "expected number seroconverted");
  Var(V,  6, "cumulative field infected");
  Var(V,  7, "cumulative housed infected");
  Var(V,  8, "cumulative poor condition infected");
  Var(V,  9, "cumulative total infected");
  Var(V, 10, "infections");
  for (i = 0; i < INTRO_MONTHS+1; i++)
    Var(V, V0+i, "P_seropos"); 
}

void Parameters(int trt, int ind, double *p, struct parameterset *Q) 
{
  Par(Q,  0, Real, "beta_{field}", Exponential, 0.02);
  Par(Q,  1, Real, "beta_{housed}", Gamma, 2.0, 0.1);
  Par(Q,  2, Real, "a", Uniform, 1.e-3, 100.);
  Par(Q,  3, Real, "k1", Exponential, 24.0);
  Par(Q,  4, Real, "k2", Exponential, 24.0);
  Par(Q,  5, Real, "L", Uniform, 0., 60.);
  Par(Q,  7, Int,  "tau", Uniform, 1, 4);

  Par(Q, 10, Real, "beta_{housed}/beta_{field}", Der);
  Par(Q, 12, Real, "E[inf_{field}]", Der);
  Par(Q, 13, Real, "E[inf_{housed}]", Der);
  Par(Q, 14, Real, "E[inf]", Der);
  Par(Q, 16, Real, "s_{median,1}", Der);
  Par(Q, 17, Real, "s_{median,2}", Der);
  Par(Q, 18, Real, "s_{95,1}", Der);
  Par(Q, 19, Real, "s_{95,2}", Der);
}

int Model(const uint trt, const uint ind, const uint nmonths, double *p, double **V, TimePoints *Data) 
{
  int i, j, n, m;
  double *Einf_field, *Einf_housed;
  double **p_seropos = V, pinf;
 
  double beta_field = p[0];
  double beta_housed = p[1];
  double S_a = p[2];
  double S_b1 = S_a/p[3];
  double S_b2 = S_a/p[4];
  double L = p[5];
  int tau = lrint(p[7]);

  double *num_infectious = doublevector(nmonths);
  double *num_infected = doublevector(nmonths);
  // double *num_seroconverted = doublevector(nmonths);
  double *L_cdf = doublevector(nmonths);
  // double *S_cdf = doublevector(nmonths);

  parameters_t pars;
  pars.Lambda = doublevector(nmonths);
  pars.lambda = doublevector(nmonths);
  pars.B = doublevector(nmonths);
  pars.b = doublevector(nmonths);
  pars.W4 = doublematrix(nmonths, nmonths);
  pars.W2 = doublevector(nmonths);
  pars.a = S_a;

  // cumulative distribution of latent period
  for (m = 0; m < nmonths; m++) {
    // L_cdf[m] = gsl_cdf_gamma_P((double)m, L_a, 1./L_b);
    L_cdf[m] = (m < L? 0: 1);
    // S_cdf[m] = gsl_cdf_gamma_P((double)m, pars.a, 1./pars.b);
  }

  // calculate expected number of infectious ewes and force of infection, lambda
  // the two donor ewes and the lamb are infected at start of experiment
  int donors_lamb[] = {0, 1, 44};
  for (i = 0; i < NUMBER_DONORS+1; i++) {
    j = donors_lamb[i];
    for (m = sheep_intro_month[j]; m < sheep_last_month[j]; m++) {
      num_infected[m]++;
      if (j == 0)
        num_infectious[m]++;
      else
        num_infectious[m] += L_cdf[m-sheep_intro_month[j]];
    }
  }

  for (m = 0; m < nmonths; m++) {
    // calculate force of infection
    // ewes are housed for tau months (1 to 4) ending in March
    if (4-tau <= (m+4)%12 && (m+4)%12 <= 3)
      // transmission rate for housed ewes in winter
      pars.lambda[m] = num_infectious[m]*beta_housed;
    else
      // field transmission rate
      pars.lambda[m] = num_infectious[m]*beta_field;

    // cumulative lambda used for prob of uninfected
    if (m > 1)
        pars.Lambda[m] = pars.Lambda[m-1]+pars.lambda[m-1];

    for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
      if (use_sheep[j] && negatives[j] > ind && sheep_intro_month[j] <= m && m < sheep_last_month[j]) {
        // probability this sheep is infected during this month m
        pinf = (1.-exp(-pars.lambda[m]))*u(sheep_intro_month[j], m, pars.Lambda);
        for (n = m+1; n < sheep_last_month[j]; n++) {
            num_infected[n] += pinf;
            num_infectious[n] += pinf*L_cdf[n-1-m];
            // num_seroconverted[i] += pinf*S_cdf[n-1-m];
        }
      }
  }

  // monthly seroconversion rates
  for (m = 0; m < nmonths; m++)
    if (m == 22)
      // increased seroconversion rate in winter 1982 due to poor condition and 
      pars.b[m] = S_b2;
    else
      pars.b[m] = S_b1;
  // cumulative seroconversion rate
  for (m = 1; m < nmonths; m++)
    pars.B[m] = pars.B[m-1] + pars.b[m-1];
  
  // W4 and W2 per month
  for (m = 0; m < nmonths; m++) {
    for (n = 0; n < m; n++)
      pars.W4[n][m] = w(n, m, n+1, m+1, &pars) 
                    - w(n, m, n+1, m,   &pars) 
                    - w(n, m, n,   m+1, &pars)  
                    + w(n, m, n,   m,   &pars);
    pars.W2[m] = w(m, m, m, m+1, &pars);
  }
 
  // calculate the probability that a randomly selected sheep seroconverts in month m
  for (j = 0; j < INTRO_MONTHS; j++)
    for (m = intro_months[j]; m < nmonths; m++)
      p_seropos[m][V0+j] = Prob_seroconvert(intro_months[j], m, &pars);

  // 3 sheep were born in month 12, were removed in month 20 and retested in month 36, all were negative
  j = 13;
  // save lambda
  double *lambda_s = doublevector(nmonths);
  double *Lambda_s = doublevector(nmonths);
  for (m = sheep_last_month[j]; m < sheep_last_test_month[j]; m++) {
    lambda_s[m] = pars.lambda[m];
    Lambda_s[m] = pars.Lambda[m];
  }

  // set lambda=0 for the months from removal (last_month) until testing and recalculate Lambda
  for (m = sheep_last_month[j]; m < sheep_last_test_month[j]; m++) 
    pars.lambda[m] = 0;
  for (m = sheep_last_month[j]+1; m < sheep_last_test_month[j]; m++)
    pars.Lambda[m] = pars.Lambda[m-1];

  // recalculate prob of seroconverting for these 3 sheep
  for (m = sheep_last_month[j]; m < sheep_last_test_month[j]; m++) {
    for (n = sheep_last_month[j]; n < m; n++)
      pars.W4[n][m] = w(n, m, n+1, m+1, &pars) 
                    - w(n, m, n+1, m,   &pars) 
                    - w(n, m, n,   m+1, &pars)  
                    + w(n, m, n,   m,   &pars);
    pars.W2[m] = w(m, m, m, m+1, &pars);
  }
  for (m = sheep_intro_month[j]; m < sheep_last_test_month[j]; m++)
    p_seropos[m][V0+sheep_var[j]] = Prob_seroconvert(sheep_intro_month[j], m, &pars);

  // restore lambda
  for (m = sheep_last_month[j]; m < sheep_last_test_month[j]; m++) {
    pars.lambda[m] = lambda_s[m];
    pars.Lambda[m] = Lambda_s[m];
  }
  free(lambda_s);
  free(Lambda_s);

  /******************** OUTPUT ***********************/
  if (Data->mode == OUTPUT) {
    for (m = 0; m < nmonths; m++)
      V[m][0] = m;

    // expected number of infected sheep in month m
    for (m = 0; m < nmonths-1; m++)
      V[m][1] = num_infected[m];
    V[nmonths-1][1] = V[nmonths-2][1];

    // expected number of infectious sheep in month m
    for (m = 0; m < nmonths-1; m++)
      V[m][2] = num_infectious[m];
    V[nmonths-1][2] = V[nmonths-2][2];

    // expected number of seroconverted sheep in month m
    // for (m = 0; m < nmonths-1; m++)
    //   V[m][5] = num_seroconverted[m];
    // V[nmonths-1][5] = V[nmonths-2][5];

    // expected number of seroconversions by the end of each month m
    // (which is the start of month m+1)
    // done this way because earlier observed seroconversions are recorded at start of month
    V[0][4] = 0;
    for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
      if (use_sheep[j] && negatives[j] > ind)
        for (m = sheep_intro_month[j]; m < sheep_last_month[j]; m++)
          V[m+1][4] += p_seropos[m][V0+sheep_var[j]];

    // cumulative expected number of seroconversions by end of month m
    // (which is the start of month m+1)
    V[0][3] = 0;
    for (m = 0; m <nmonths-1; m++)
      V[m+1][3] = V[m][3]+V[m+1][4];

    // expected number of infection events in each month m partitioned by route
    Einf_field = doublevector(nmonths);
    Einf_housed = doublevector(nmonths);
    for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
      if (use_sheep[j] && negatives[j] > ind) {
        for (m = sheep_intro_month[j]; m < sheep_last_month[j]; m++)
          if (4-tau <= (m+4)%12 && (m+4)%12 <= 3)
            Einf_housed[m] += (1.-exp(-pars.lambda[m]))*u(sheep_intro_month[j], m , pars.Lambda);
          else
            Einf_field[m] += (1.-exp(-pars.lambda[m]))*u(sheep_intro_month[j], m, pars.Lambda);
      }
    for (m = 0; m < nmonths; m++) 
      V[m][10] = Einf_field[m]+Einf_housed[m];

    for (m = 1; m < nmonths; m++) {
      Einf_field[m] += Einf_field[m-1];
      Einf_housed[m] += Einf_housed[m-1];
    }
    for (m = 0; m < nmonths; m++) {
      V[m][6] = Einf_field[m];
      V[m][7] = Einf_housed[m];
      V[m][9] = V[m][6]+V[m][7];
    }
    free(Einf_field);
    free(Einf_housed);
  }

  free(L_cdf);
  // free(S_cdf);
  free(num_infected);
  free(num_infectious);
  // free(num_seroconverted);
  free(pars.lambda);
  free(pars.Lambda);
  free(pars.b);
  free(pars.B);
  free(pars.W4[0]); 
  free(pars.W4);
  free(pars.W2);
  return SUCCESS;
}

double P_pos(int j, int s1, int s2, double **V)
{
  int m;
  double sum = 0;
  for (m = s1; m < s2; m++)
    sum += V[m][V0+sheep_var[j]];
  return sum;
}

double logLikelihood(const uint trt, const uint ind, const double *p, double **V, const TimePoints *Data)
{
  int j, pos_idx;
  double lnL = 0;

  for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
    if (use_sheep[j] && negatives[j] > ind) {
      pos_idx = sheep_pos_index[j];
      if (pos_idx > 0)
        // ewe is first seropositive at start of month m, therefore it seroconverted
        // sometime in the months from the last test to the month prior to seroconverting
        lnL += log(P_pos(j, lrint(Data->t[pos_idx-1][0]), sheep_pos_month[j], V));
      else
        // ewe never seroconverted
        lnL += log(1. - P_pos(j, sheep_intro_month[j], sheep_last_test_month[j], V));
    }
  return lnL;
}

void OutputModel(int trt, int ind, double *Output, double *V)
{
  int i;
  for (i = 1; i < V0+INTRO_MONTHS+1; i++)
    Output[i] = V[i];
}

void OutputData(int trt, int ind, double *Output, double *Var, double *Data, uint *value, int index)
{
  Output[3] = Data[47];
  Output[4] = Data[48];
  Output[5] = Data[49];
  // remove dam to lamb infection from data to match model output
  if (lrint(Data[0]) == 3)
    Output[4] = 0;
  if (lrint(Data[0]) >= 3)
    Output[3] -= 1;
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
  F->n_func = 3;
  F->function_list = (function_ptr*) malloc(F->n_func*sizeof(function_ptr));
  F->function_list[0] = f0;
  F->function_list[1] = f1;
  F->function_list[2] = f2;
  // F->function_list[3] = f3;
}

double u(int m0, int m, double *Lambda)
{
  // The probability a sheep is uninfected at the start of month m 
  // given it was introduced into the flock in month m0
  if (m == m0)
    return 1;
  else
    return exp(-Lambda[m]+Lambda[m0]);
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

double w(int n, int m, int t, int s, parameters_t *p)
{
  double a = p->a;
  double l = p->lambda[n];
  double b = p->b[n];
  double r = (b-l)/b;
  double Z, z, v;

  if (m-n <= 1 && t == m && s == m)
    return 0;

  Z = p->B[m] - p->B[n] + b*(n-t) + p->b[m]*(s-m);
  if (Z == 0)
    return 0;

  z = r*Z;
  if (z <= 0)
    v = exp(a*log(Z) - l*Z/b + ln_inc_gamma_star(a, z));
  else
    v = exp(-a*log(r) - l*Z/b + log(P_inc_gamma(a, z)));

  return exp(-l*(t-n))*(v - P_inc_gamma(a, Z));
}

double Prob_seroconvert(int m0, int m, parameters_t *p)
{
  int n;
  double P_m = 0;

  for (n = m0; n < m; n++)
      P_m += u(m0, n, p->Lambda)*p->W4[n][m];
  P_m -= u(m0, m, p->Lambda)*p->W2[m];

  return min(max(0, P_m), 1);
}

void Residual(int trt, int ind, double *R, double *V, double *p, double *Data, uint *value) 
{
  R[3] = V[3]-Data[47];
  R[4] = V[4]-Data[48];
}

void DerivedParameters(int trt, int ind, int nmonths, double *p, double **V, TimePoints *Data) 
{
  // int i, j, m;
  // double pinf;
  double beta_field = p[0];
  double beta_housed = p[1];
  double S_a = p[2];
  double S_b1 = S_a/p[3];
  double S_b2 = S_a/p[4];
  // double L_mu = p[5];
  // double L_sd = p[6];
  // double *num_infectious = doublevector(nmonths);
  // double *Lambda = doublevector(nmonths);
  // double *lambda = doublevector(nmonths);
  // double *B = doublevector(nmonths);
  // double *b = doublevector(nmonths);
  // double *L_cdf = doublevector(nmonths);

  p[10] = beta_housed/beta_field;
  p[16] = P_inv_gamma(0.5, S_a, 1./S_b1);
  p[17] = P_inv_gamma(0.5, S_a, 1./S_b2);
  p[18] = P_inv_gamma(0.95, S_a, 1./S_b1);
  p[19] = P_inv_gamma(0.95, S_a, 1./S_b2);
  return;

  // p[6] = gsl_cdf_exponential_Pinv(0.5, mu);
  // p[22] = gsl_cdf_exponential_Pinv(0.5, mu2);
  // p[15] = gsl_cdf_exponential_Pinv(0.025, mu);
  // p[16] = gsl_cdf_exponential_Pinv(0.975, mu);
  // p[20] = 6./16.*(1.-exp(-(0.75*beta_field)*(fmax(0, 4*12-L_mu))));
  // p[23] = 6./16.*(1.-exp(-(0.25*beta_housed)*(fmax(0, 4*12-L_mu))));

  // for (m = 0; m < nmonths; m++)
  //    L_cdf[m] = (m < L_mu? 0: 1);

  // // calculate expected number of infectious ewes and force of infection, lambda
  // // the two donor sheep are infectious at start of experiment
  // for (j = 0; j < 2; j++)
  //   for (m = 0; m < sheep_last_month[j]; m++)
  //     num_infectious[m]++;
  // // ewe 44 - infected by dam
  // for (m = sheep_intro_month[44]; m < sheep_last_month[44]; m++)
  //   num_infectious[m] += L_cdf[m-sheep_intro_month[44]];

  // for (m = 0; m < nmonths; m++) {
  //   // calculate force of infection
  //   if ((m+3)%12 < 3)
  //     // transmission rate for housed ewes (in winter Jan-Mar)
  //     lambda[m] = num_infectious[m]*beta_housed;
  //   else
  //     // field transmission rate
  //     lambda[m] = num_infectious[m]*beta_field;

  //   // cumulative lambda for prob of uninfected
  //   if (m > 1)
  //       Lambda[m] = Lambda[m-1]+lambda[m-1];

  //   for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
  //     if (use_sheep[j] && sheep_intro_month[j] <= m && m < sheep_last_month[j]) {
  //       pinf = (1.-exp(-lambda[m]))*u(sheep_intro_month[j], m, Lambda);
  //       for (i = m; i < sheep_last_month[j]; i++)
  //           num_infectious[i] += pinf*L_cdf[i-m];
  //     }
  // }

  // p[12] = 0;
  // p[13] = 0;
  // for (j = NUMBER_DONORS; j < NUMBER_SHEEP; j++)
  //   if (use_sheep[j]) {
  //     for (m = sheep_intro_month[j]; m < sheep_last_month[j]; m++)
  //       if ((m+3)%12 < 3)
  //         p[13] += (1.-exp(-lambda[m]))*u(sheep_intro_month[j], m, Lambda);
  //       else
  //         p[12] += (1.-exp(-lambda[m]))*u(sheep_intro_month[j], m, Lambda);
  //   }
  // p[14] = p[12]+p[13];
  
  // free(lambda);
  // free(Lambda);
  // free(B);
  // free(b);
  // free(num_infectious);
  // free(L_cdf);
}

void WAIC(int trt, int ind, double **lnL, double **V, double *p, const TimePoints *Data) {}
void PredictData(int trt, int ind, double *Output, double *V, double *p, gsl_rng *stream) {}
void SimulateData(int trt, int ind, double *Output, double *V, double *p, double *Data, uint *value, gsl_rng *stream) {}
double timestep(void) {return 1;}
void GlobalParameters() {}
double UserPDF(double x) {return 0;}
void HyperParameters(Treatments T, Hyperparameters *H) {}
void SaturatedModel(int trt, int ind, double **V, double *p, const TimePoints *Data) {}

