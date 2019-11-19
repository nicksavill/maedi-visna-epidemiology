#include <func.h>

void ReadParam(char *type, void *var, char *varname) 
{
  char parname[1000];
  char format[4][10] = { "%d", "%f", "%lf", "%s" };
  int typeno;
  int eofcheck = 0;
  static int first = TRUE;
  time_t t;		
  FILE *parfile, *parlog;
	
  parfile = fopen(parameterfilename, "r");
  if (parfile == NULL) FileError("parameter", parameterfilename);

  if (strcmp(type, "int") == 0)
    typeno = 0;
  else if (strcmp(type, "float") == 0)
    typeno = 1;
  else if (strcmp(type, "double") == 0)
    typeno = 2;
  else if (strcmp(type, "char*") == 0)
    typeno = 3;
  else
    Error("Can't read data type\n");
	
  while ((eofcheck = fscanf(parfile, "%s", parname)) != EOF)
    if (!strcmp(parname, varname)) break;
  if (eofcheck != EOF)
    fscanf(parfile, format[typeno], (void*) var);

  fclose(parfile);
	
  t = time(NULL);
  if (strcmp(varname, "expno") != 0 && 
       strcmp(varname, "parameterlogfilename") != 0) {
    if (first)
      parlog = fopen(parameterlogfilename, "a");
    else
      parlog = fopen(parameterlogfilename, "a");
    if (parlog == NULL) FileError("parameter log", 
				     parameterlogfilename);
	
    if (first) {
      fprintf(parlog, "%s", asctime(localtime(&t)));
      first = FALSE;
    }

    switch(typeno) {
    case 0:
      fprintf(parlog, "%s\t %d\n", varname, *(int*) var);
      break;
    case 1:
      fprintf(parlog, "%s\t %e\n", varname, *(float*) var);
      break;
    case 2:
      fprintf(parlog, "%s\t %e\n", varname, *(double*) var);
      break;
    case 3:
      fprintf(parlog, "%s\t %s\n", varname, (char*) var);
      break;
    }

    fclose(parlog);
  }
}

double **doublematrix(int sizex, int sizey)
{
  int i;
  double **m;
  //  printf("double2 %d %d\n", sizex, sizey);

  m = (double **) calloc(sizex, sizeof(double*));
  if (!m) {
    printf("%d %d\n", sizex, sizey);
    Error("Can't allocate memory in doublematrix()");
  }

  m[0] = (double *) calloc(sizex * sizey, sizeof(double));
  if (!m[0]) Error("Can't allocate memory in doublematrix()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  return m;
}

int **integermatrix(int sizex, int sizey)
{
  int i;
  int **m;
  
  m = (int**) calloc(sizex, sizeof(int*));
  if (!m) Error("Can't allocate memory in integermatrix()");
  m[0] = (int*) calloc(sizex * sizey, sizeof(int));
  if (!m[0]) Error("Can't allocate memory in integermatrix()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  return m;
}

bool **boolmatrix(int sizex, int sizey)
{
  int i;
  bool **m;
  
  m = (bool**) calloc(sizex, sizeof(bool*));
  if (!m) Error("Can't allocate memory in boolmatrix()");
  m[0] = (bool*) calloc(sizex * sizey, sizeof(bool));
  if (!m[0]) Error("Can't allocate memory in boolmatrix()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  return m;
}

uint **uintegermatrix(int sizex, int sizey)
{
  int i;
  uint **m;
	
  m = (uint**) calloc(sizex, sizeof(uint*));
  if (!m) Error("Can't allocate memory in uintegermatrix()");
  m[0] = (uint*) calloc(sizex * sizey, sizeof(uint));
  if (!m[0]) Error("Can't allocate memory in uintegermatrix()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  return m;
}

char **charmatrix(int sizex, int sizey)
{
  int i;
  char **m;
	
  m = (char**) malloc(sizex * sizeof(char*));
  if (!m) Error("Can't allocate memory in charmatrix()");
  m[0] = (char*) malloc(sizex * sizey * sizeof(char));
  if (!m[0]) Error("Can't allocate memory in charmatrix()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  return m;
}

double ***doublematrix3(int sizex, int sizey, int sizez)
{
  int i, j;
  double ***m;
  double *lastaddress;

  m = (double ***) calloc(sizex, sizeof(double**));
  if (!m) Error("Can't allocate memory in doublematrix3()");

  m[0] = (double **) calloc(sizex * sizey, sizeof(double*));
  if (!m[0]) Error("Can't allocate memory in doublematrix3()");
  for (i = 1; i < sizex; i++)
    m[i] = m[i-1] + sizey;

  m[0][0] = (double *) calloc(sizex * sizey * sizez, sizeof(double));
  if (!m[0][0]) Error("Can't allocate memory in doublematrix3()");
  for (j = 1; j < sizey; j++)
    m[0][j] = m[0][j-1] + sizez;

  lastaddress = m[0][sizey-1];
  for (i = 1; i < sizex; i++)
    for (j = 0; j < sizey; j++) {
      m[i][j] = lastaddress + sizez;
      lastaddress = m[i][j];
    }

  return m;
}

bool *boolvector(int size)
{
  bool *m;
  //  printf("bool %d\n", size);

  m = (bool*) calloc(size,  sizeof(bool));
  if (!m) {
    fprintf(stderr, "Can't allocate %d memory in boolvector()", size);
    Error("");
  }

  return m;
}

double *doublevector(int size)
{
  double *m;
  //  printf("double %d\n", size);

  m = (double*) calloc(size,  sizeof(double));
  if (!m) {
    fprintf(stderr, "Can't allocate %d memory in doublevector()", size);
    Error("");
  }

  return m;
}

double complex *doublecomplexvector(int size)
{
  double complex *m;
  //  printf("double %d\n", size);

  m = (double complex*) calloc(size,  sizeof(double complex));
  if (!m) {
    fprintf(stderr, "Can't allocate %d memory in doublecomplexvector()", size);
    Error("");
  }

  return m;
}

int *integervector(int size)
{
  int *m;
  
  m = (int*) calloc(size,  sizeof(int));
  if (!m) Error("Can't allocate memory in integervector()");

  return m;
}

uint *uintegervector(int size)
{
  uint *m;
  
  m = (uint*) calloc(size,  sizeof(uint));
  if (!m) Error("Can't allocate memory in uintegervector()");

  return m;
}

char *charvector(int size)
{
  char *m;
	
  m = (char*) calloc(size,  sizeof(char));
  if (!m) Error("Can't allocate memory in integervector()");

  return m;
}

void Error(char* message)
{
  fprintf(stderr, "%s\n", message);
  exit(1);
}

void FileError(char *type, char *filename)
{
  char message[50];

  sprintf(message, "Can't open %s file: %s", type, filename);
  Error(message);
}


void Int_HeapSort( int *x, int n )
{
  int i, j, k, l;
  int tmpx;

  if ( n < 1 ) return;

  l = n / 2 + 1;
  k = n;

  for (;;) {
    if ( l > 1 ) {
      l--;
      tmpx = x[l];
    }
    else {
      tmpx = x[k];
      x[k] = x[1];
      k--;
      if ( k <= 1 ) {
	x[1] = tmpx;
	for ( i = 1; i <= n/2; i++ ) {
	  tmpx = x[i];
	  x[i] = x[n-i+1];
	  x[n-i+1] = tmpx;
	}
	return;
      }
    }
    i = l;
    j = l+l;
    
    while ( j <= k ) {
      if ( j < k && x[j] > x[j+1] ) j++;
      if ( tmpx > x[j] ) {
	x[i] = x[j];
	i = j;
	j += j;
      }
      else
	j = k + 1;
    }
    x[i] = tmpx;
  }
}

void HeapSort(double *x, int n)
{
  int i, j, k, l;
  double tmpx;

  if (n < 1) return;

  l = n / 2 + 1;
  k = n;

  for (;;) {
    if (l > 1) {
      l--;
      tmpx = x[l];
    }
    else {
      tmpx = x[k];
      x[k] = x[1];
      k--;
      if (k <= 1) {
	x[1] = tmpx;
	for (i = 1; i <= n/2; i++) {
	  tmpx = x[i];
	  x[i] = x[n-i+1];
	  x[n-i+1] = tmpx;
	}
	return;
      }
    }
    i = l;
    j = l+l;
    
    while (j <= k) {
      if (j < k && x[j] > x[j+1]) j++;
      if (tmpx > x[j]) {
	x[i] = x[j];
	i = j;
	j += j;
      }
      else
	j = k + 1;
    }
    x[i] = tmpx;
  }
}

void circmeanvar(double *x, int n, double *mean, double *var)
{
  int i;
  complex m1 = 0+I*0;

  for (i = 0; i < n; i++)
    m1 += cos(x[i])+I*sin(x[i]);
  m1 /= n;

  /* *mean = fmod(2.*M_PI+carg(m1), 2.*M_PI)/(2.*M_PI)*24.; */
  /* //  *var = -2*log(cabs(m1))/(4.*M_PI*M_PI)*24.*24.; */
  /* *var = (1.-cabs(m1))/(4.*M_PI*M_PI)*24.*24.; */

  *mean = fmod(2.*M_PI+carg(m1), 2.*M_PI);
  //  *var = -2*log(cabs(m1));
  *var = 1.-cabs(m1);
}

void meanvar(double *x, int n, double *m, double *s2)
{
  int i, c = 0;
  double a, b = 0;

  *m = 0;
  *s2 = 0;
  if (n == 0) return;
  if (n == 1) {
    *m = x[0];
    return;
  }
  for (i = 0; i < n; i++)
    if (isfinite(x[i])) {
      *m += x[i];
      c++;
    }
  *m /= (double)c;

  for (i = 0; i < n; i++)
    if (isfinite(x[i])) {
      a = (x[i]-(*m));
      b += a;
      *s2 += a*a;
    }
  *s2 = ((*s2)-b*b/(double)c)/(double)(c-1);
}

void meanvarl(long double *x, int n, long double *m, long double *s2)
{
  int i;
  long double a, b = 0;

  *m = 0;
  *s2 = 0;
  if (n == 0) return;
  if (n == 1) {
    *m = x[0];
    return;
  }
  for (i = 0; i < n; i++)
    *m += x[i];
  *m /= (long double)n;

  for (i = 0; i < n; i++) {
    a = (x[i]-(*m));
    b += a;
    *s2 += a*a;
  }
  *s2 = ((*s2)-b*b/(long double)n)/(long double)(n-1);
}

void wmeanvar(long double *x, int n, long double *m, long double *ss2)
{
  // mean and variance for very large numbers
  // the variance is scaled by 1/mean^2
  // on return you need to multiply ss2 by mean^2 to get the variance

  int i;
  long double a, b = 0;

  *m = 0;
  *ss2 = 0;
  if (n == 0) return;
  if (n == 1) {
    *m = x[0];
    return;
  }
  for (i = 0; i < n; i++)
    *m += x[i];

  for (i = 0; i < n; i++) {
    a = x[i]/(*m)-1.;
    b += a;
    *ss2 += a*a;
  }
  *ss2 = ((*ss2)-b*b/(long double)n)/(long double)(n-1);
}

double mean(double *x, int s, int n)
{
  int i;
  double m = 0;

  if (n - s <= 0) {
    fprintf(stderr, "s=%d n=%d\n", s, n);
    Error("range wrong for mean");
  }
  for (i = s; i < n; i++)
    m += x[i];
  return m/(double)(n-s);
}

double poissonCMF(const double x, const double rate, int type)
{
  if (rate <= 0) Error("rate <=0 in poissonCDF");
  if (x < 0)
    return -INFINITY;
  else if (type == ABOVE_DETECTION)
    return log(gsl_cdf_poisson_Q((uint)x, rate));
  else if (type == BELOW_DETECTION)
    return log(gsl_cdf_poisson_P((uint)x, rate));
  else
    Error("use ABOVE_DETECTION or BELOW_DETECTION for cumulative poisson"
	  "distribution");
  return -INFINITY;
}

double gamma_scalePDF(const double x, const double shape, const double scale)
{
  if (shape <= 0) Error("shape <=0 in gamma_scalePDF");
  if (scale <= 0) Error("scale <=0 in gamma_scalePDF");
  if (x < 0)
    return -INFINITY;
  else
    return (shape-1.)*log(x) - x/scale;
}

double betaPDF(const double x, const double alpha, const double beta)
{
  if (alpha <= 0) Error("alpha <=0 in betaPDF");
  if (beta <= 0) Error("beta <=0 in betaPDF");
  if (x < 0 || x > 1)
    return -INFINITY;
  else
    return (alpha-1.)*log(x)+(beta-1.)*log1p(-x);
}

double exponentialPDF(const double x, const double mean)
{
  if (mean < 0) Error("mean < 0 in exponentialPDF");
  if (x < 0)
    return -INFINITY;
  else
    return -x/mean;
}

double poissonPMF(const double x, const double mean)
{
  if (mean < 0) Error("mean < 0 in poissonPMF");
  if (x < 0) Error("x < 0 in poissonPMF");
  if (mean == 0) 
    if (x == 0)
      return 0;
    else
      return -INFINITY;
  else if (x == 0)
    return -mean;
  else 
    return x*log(mean)-gsl_sf_lnfact((uint)x);
}

double geometricPDF(const double x, const double p)
{
  if (p < 0) Error("p < 0 in geometricPDF");
  if (p > 1) Error("p > 1 in geometricPDF");
  if (x < 0)
    return -INFINITY;
  else
    return x*log1p(-p);
}

double trnormalPDF(const double x, const double mean, const double var)
{
  if (var < 0) Error("var < 0 in trnormalPDF");
  if (x < 0)
    return -INFINITY;
  else {
    double u = x-mean;
    return -u*u/(2.*var);
  }
}

double normalCDF(const double x, const double mean, const double var, int type)
{
  if (var < 0) Error("var < 0 in normalCDF");
  if (type == ABOVE_DETECTION)
    return log(gsl_cdf_gaussian_Q(x-mean, sqrt(var)));
  else if (type == BELOW_DETECTION)
    return log(gsl_cdf_gaussian_P(x-mean, sqrt(var)));
  else
    Error("use ABOVE_DETECTION or BELOW_DETECTION for cumulative normal"
	  "distribution");
  return -INFINITY;
}

double normalPDF(const double x, const double mean, const double var)
{
  if (var < 0) Error("var < 0 in normalPDF");
  double u = x-mean;
  return -u*u/(2.*var);
}

double binomialPMF(const double k, const double n, const double p)
{
  /* 
     k  no. success
     n  no. trials
     p  prob success
  */
  
  if (isnan(p)) Error("p is NaN in binomialPMF");
  if (n < 0) Error("n < 0 in binomialPMF");
  if (p < 0) Error("p < 0 in binomialPMF");
  if (p > 1) Error("p > 1 in binomialPMF");
  if (k < 0 || n < k) return -INFINITY;
  if (p == 0)
    return (k == 0)? 0: -INFINITY;
  else if (p == 1)
    return (k == n)? 0: -INFINITY;
  else
    return k*log(p)+(n-k)*log1p(-p);
}

double multinomialPMF(const unsigned int K, const double *r, const double *p)
{
  /* K categories
     r[0..K-1] successes
     p[0..K-1] probabilites
  */
  unsigned int i;
  double l = 0;
  
  for (i = 0; i < K; i++) {
    if (isnan(p[i])) Error("p is NaN in multinomialPMF");
    if (r[i] < 0) Error("r < 0 in multinomialPMF");
    if (p[i] < 0) Error("p < 0 in multinomialPMF");
    if (p[i] > 1) Error("p > 1 in multinomialPMF");
    l += p[i];
  }
  if (lrint(l*1e6) > 1000000) {
    printf("sum = %.20f\n", l);
    Error("total p > 1 in multinomialPMF");
  }
  
  l = 0;
  for (i = 0; i < K; i++)
    if (p[i] == 0) {
      if (r[i] != 0) return -INFINITY;
    }
    else if (r[i] != 0)
      l += r[i]*log(p[i]);
  return l;
}

double uniformPDF(const double x, const double a, const double b)
{
  if (b < a) Error("b < a in uniformPDF");
  if (x < a || x > b)
    return -INFINITY;
  else
    return -log(b-a);
}

double logisticPDF(const double x, const double mean, const double var)
{
  if (var < 0) Error("var < 0 in logisticPDF");
  const double u = (x-mean)/sqrt(var);
  return -u-2.*log(1+exp(-u));
}

double betabinomialPMF(const double x, const double n, const double alpha,
		       const double beta)
{
  if (alpha < 0) Error("alpha negative in betabinomialPMF");
  if (beta < 0) Error("beta negative in betabinomialPMF");
  if (n < 0) Error("n negative in in betabinomialPMF");
  if (x < 0) Error("x negative in in betabinomialPMF");
  if (x > n)  Error("x>n in in betabinomialPMF");
  return gsl_sf_lnbeta(n+alpha, n-x+beta) - gsl_sf_lnbeta(alpha, beta);
}

/* double MVNPDF(Matrix *x, Matrix *mu, Matrix *prec, const double logdetprec, */
/* 	      int d) */
/* { */
/*   double l = 0; */
/*   Matrix xt, z, mus, precs; */

/*   Matrix_Init(&xt, "xt", 1, d); */
/*   Matrix_Init(&z, "z", 1, 1); */
/*   Matrix_Init(&mus, "mus", d, 1); */
/*   Matrix_Init(&precs, "precs", d, d); */

/*   Matrix_Sub_Copy(&mus, mu, d); */
/*   Matrix_Sub_Copy(&precs, prec, d); */

/*   Matrix_Subtract(x, x, &mus); */
/*   Matrix_Transpose(&xt, x); */
/*   Matrix_Multiply(&xt, &xt, &precs); */
/*   Matrix_Multiply(&z, &xt, x); */

/*   l = -0.5*((double)d*LN2PI - logdetprec + z.a[0]); */
/*   Matrix_Destroy(&xt); */
/*   Matrix_Destroy(&z); */
/*   Matrix_Destroy(&mus); */
/*   Matrix_Destroy(&precs); */
/*   return l; */
/* } */

#define JMAX2 40
double rtbis(double (*func)(double, double*), double x1, double x2,
 	      double xacc, int *flag, double *p)
{
  int j;
  double dx, f, fmid, xmid, rtb;
  *flag = TRUE;
  
  f = (*func)(x1, p);
  fmid = (*func)(x2, p);
  if (isnan(f) || isnan(fmid)) Error("NaN in Bisect: widen bounds");
  if (f*fmid >= 0.0) {
    *flag = FALSE;
    printf("not bracketed\n");
    return FALSE;
  }
  rtb = f < 0.0 ? (dx = x2-x1,x1) : (dx = x1-x2,x2);
  
  for (j = 1;j<=JMAX2;j++) {
    fmid = (*func)(xmid = rtb+(dx *= 0.5), p);
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  *flag = FALSE;
  printf("Too many bisections in rtbis\n");
  return FALSE;
 }

double normalRNG(double mu, double var, const gsl_rng *stream)
{
  if (var < 0) {
    fprintf(stderr, "var<0 in normalRNG\n");
    Error("");
  }
  return mu + gsl_ran_gaussian(stream, sqrt(var));
}

uint negative_binomialRNG(double n, double p, const gsl_rng *stream)
{
  // return number of failures before n successes in independent trials with
  // probability p of success
  
  if (isnan(p)) Error("p is NaN");
  if (n < 0) {
    fprintf(stderr, "n = %f in negative_binomialRNG\n", n);
    Error("");
  }
  if (p < 0. || p > 1.) {
    fprintf(stderr, "p = %e in negative_binomialRNG\n", p);
    Error("");
  }
  // needed because GSL library enters infinite loop if np is very small
  // mean of negative binoimal is np/(1-p)
  if (n*p < 1e-3) return 0; 

  return gsl_ran_negative_binomial(stream, p, n);
}

uint binomialRNG(uint n, double p, const gsl_rng *stream)
{
  if (isnan(p)) Error("p is NaN");
  if (n < 0) {
    fprintf(stderr, "n = %d in binomialRNG\n", n);
    Error("");
  }
  if (p < 0. || p > 1.) {
    fprintf(stderr, "p = %e in binomialRNG\n", p);
    Error("");
  }

  return gsl_ran_binomial(stream, p, n);
}

void multinomialRNG(uint n, const double *p, size_t K, unsigned int *r, 
		    const gsl_rng *stream)
{
  /* K categories
     n trials
     p[0..K-1] probabilites
     r[0..K-1] successes */

  gsl_ran_multinomial(stream, K, n, p, r);
}

double uniformRNG(double a, double b, const gsl_rng *stream)
{
  if (a >= b) {
    fprintf(stderr,"Uniform: %e %e\n", a, b);
    Error("");
  }
  return gsl_ran_flat(stream, a, b);
}

int uniform_intRNG(int a, int b, const gsl_rng *stream)
{
  if (a >= b) {
    fprintf(stderr,"Uniform_int: %d %d\n", a, b);
    Error("");
  }
  return a+(int)gsl_rng_uniform_int(stream, (long int)(b-a));
}

int poissonRNG(double lambda, const gsl_rng *stream)
{
  if (lambda < 0) Error("lambda < 0 in poissonRNG");
  if (isnan(lambda)) Error("lambda is NaN in poissonRNG");
  return gsl_ran_poisson(stream, lambda);
}

void dirichletRNG(const double *alpha, size_t K, double *r,
		  const gsl_rng *stream)
{
  gsl_ran_dirichlet(stream, K, alpha, r);
}

double gammaRNG(double shape, double scale, const gsl_rng *stream)
{
  // scale is the scale parameter (theta) NOT the inverse scale (beta)
  if (shape <= 0) {
    fprintf(stderr, "gammaRNG: shape = %e should not be less than zero\n", shape);
    Error("");
  }
  if (scale <= 0) {
    fprintf(stderr, "gammaRNG: scale = %e should not be less than zero\n", scale);
    Error("");
  }
  return gsl_ran_gamma(stream, shape, scale);
} 

double chisqRNG(double nu, const gsl_rng *stream)
{
  return gsl_ran_chisq(stream, nu);
}

double scaleinvchisqRNG(double nu, double s2, const gsl_rng *stream) 
{
  return nu*s2/chisqRNG(nu, stream);
}

double invgammaRNG(double shape, double scale, const gsl_rng *stream) 
{
  return 1./gammaRNG(shape, 1./scale, stream);
}

double paretoRNG(double scale, double order, const gsl_rng *stream) 
{
  return gsl_ran_pareto(stream, order, scale);
}

double betaRNG(double alpha, double beta, const gsl_rng *stream) 
{
  if (alpha <= 0) Error("alpha<=0 in betaRNG");
  if (isinf(alpha)) Error("alpha=infinity in betaRNG");
  if (beta <= 0) Error("beta<=0 in betaRNG");
  if (isinf(beta)) Error("beta=infinity in betaRNG");
  return gsl_ran_beta(stream, alpha, beta);
}

void dirichlet_multinomialRNG(uint n, const double *alpha, size_t K, uint *r,
			      const gsl_rng *stream)
{
  double p[K];
  
  dirichletRNG(alpha, K, p, stream);  
  multinomialRNG(n, p, K, r, stream);
}

uint beta_binomialRNG(uint n, double alpha, double beta, const gsl_rng *stream)
{
  if (alpha <= 0) {
    fprintf(stderr, "beta_binomialRNG: alpha = %e should not be less than zero\n", alpha);
    Error("");
  }
  if (beta <= 0) {
    fprintf(stderr, "beta_binomialRNG: beta = %e should not be less than zero\n", beta);
    Error("");
  }
  return binomialRNG(n, betaRNG(alpha, beta, stream), stream);
}

uint beta_negative_binomialRNG(uint n, double alpha, double beta,
			       const gsl_rng *stream)
{
  if (alpha <= 0) {
    fprintf(stderr, "beta_negative_binomialRNG: alpha = %e should not be less than zero\n", alpha);
    Error("");
  }
  if (beta <= 0) {
    fprintf(stderr, "beta_negative_binomialRNG: beta = %e should not be less than zero\n", beta);
    Error("");
  }
  double r = betaRNG(alpha, beta, stream);
  //  printf("%f\n", r);
  return negative_binomialRNG((double)n, r, stream);
}

/* Matrix *wishartRNG(int d, int n, Matrix *S, Matrix *R, gsl_rng *stream) */
/* { */
/*   int i, j; */
/*   Matrix A; */
  
/*   if (n <= d-1) { */
/*     fprintf(stderr, "n=%d, d=%d\n", n, d); */
/*     Error("wishart"); */
/*   } */
  
/*   Matrix_Init(&A, "wishart", d, d); */
  
/*   // R and matlab implememtation */
/*   // A is lower triangular */
/*   for (i = 0; i < d; i++) { */
/*     A.a[i*d+i] = sqrt(chisqRNG((double)(n-i), stream)); */
/*     for (j = 0; j < i; j++) */
/*       A.a[i*d+j] = normalRNG(0., 1., stream); */
/*   } */
  
/*   // S=chol(S), R=S'A'AS=(AS)'AS */
/*   if (Matrix_Cholesky(S, S) == NULL) { */
/*     Matrix_Destroy(&A); */
/*     return NULL; */
/*   } */
/*   Matrix_Multiply(S, &A, S); */
/*   Matrix_Transpose(&A, S); */
/*   Matrix_Multiply(R, &A, S); */

/*   Matrix_Destroy(&A); */
/*   return R; */
/* } */

/* void MVNRNG(int d, Matrix *mu, Matrix *cholS, Matrix *r, gsl_rng *stream)  */
/* { */
/*   vdRngGaussianMV(VSL_METHOD_DGAUSSIANMV_BOXMULLER2, stream, 1, r->a, d,  */
/* 		  VSL_MATRIX_STORAGE_FULL, mu->a, cholS->a); */
/* } */

Matrix *MatrixInit(Matrix *A, int n, int m, int *elements) 
{
  int i, j, k = 0;
  if (n < 1) Error("Init: > 0 rows for matrix");
  if (m < 1) Error("Init: > 0 columns for matrix");

  A->n = n;
  A->m = m;
  A->a = integermatrix(n, m);
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      A->a[i][j] = elements[k++];
  return A;
}

Matrix *MatrixDestroy(Matrix *A)
 {
   free(A->a);
   A->n = A->m = 0;
   return A;
 }

/* Matrix *Matrix_Init(Matrix *A, char *name, int n, int m)  */
/* { */
/*   if (n < 1) Error("Init: > 0 rows for matrix"); */
/*   if (m < 1) Error("Init: > 0 columns for matrix"); */

/*   strcpy(A->name, name); */

/*   A->n = n; */
/*   A->m = m; */
/*   A->a = doublevector(n*m); */
/*   return A; */
/* } */

/* Matrix *Matrix_Copy(Matrix *X, Matrix *A) */
/* { */
/*   int i; */
/*   if (X != A) { */
/*     if (X->n != A->n || X->m != A->m) Error("Copy: matrices not of same order"); */
/*     for (i = 0; i < X->n*X->m; i++) */
/*       X->a[i] = A->a[i]; */
/*   } */
/*   return X; */
/* } */

/* Matrix *Matrix_Sub_Copy(Matrix *X, Matrix *A, int n) */
/* { */
/*   int i, j; */
/*   for (i = 0; i < MIN(n, X->n); i++) */
/*     for (j = 0; j < MIN(n, X->m); j++) */
/*       X->a[i*X->m+j] = A->a[i*A->m+j]; */
/*   return X; */
/* } */

/* Matrix *Matrix_Transpose(Matrix *X, Matrix *A) */
/* { */
/*   int i, j; */
/*   double *a = doublevector(X->n*X->m); */

/*   if (X->n != A->m || X->m != A->n) Error("matrices not transposable"); */
/*   for (i = 0; i < A->n; i++) */
/*     for (j = 0; j < A->m; j++) */
/*       a[j*X->m+i] = A->a[i*A->m+j]; */

/*   for (i = 0; i < X->n*X->m; i++) */
/*     X->a[i] = a[i]; */

/*   free(a); */
/*   return X; */
/* } */

/* Matrix *Matrix_Scalar_Multiply(Matrix *X, double lambda, Matrix *A) */
/* { */
/*   int i; */
/*   for (i = 0; i < X->n*X->m; i++) */
/*     X->a[i] = lambda*A->a[i]; */
/*   return X; */
/* } */

/* Matrix *Matrix_Add(Matrix *X, Matrix *A, Matrix *B) */
/* { */
/*   int i; */
/*   if (X != A && (X->n != A->n || X->m != A->m) || */
/*       X != B && (X->n != B->n || X->m != B->m)) { */
/*     Matrix_Print(A); */
/*     Matrix_Print(B); */
/*     Error("Add: matrices not of same order"); */
/*   } */
/*   for (i = 0; i < X->n*X->m; i++) */
/*     X->a[i] = A->a[i] + B->a[i]; */
/*   return X; */
/* } */

/* Matrix *Matrix_Add_Diagonal(Matrix *A, Matrix *B) */
/* { */
/*   // Add the diagonal of A to B, return B. A must be diagonal */
/*   int i, j, n = A->n; */
/*   if (n != A->m || n != B->n) { */
/*     Matrix_Print(A); */
/*     Matrix_Print(B); */
/*     Error("Add Diagonal: matrices not square or of same order"); */
/*   } */
/*   for (i = 0; i < n; i++) { */
/*     j = i*n+i; */
/*     B->a[j] = A->a[j] + B->a[j]; */
/*   } */
/*   return B; */
/* } */

/* Matrix *Matrix_Subtract(Matrix *X, Matrix *A, Matrix *B) */
/* { */
/*   int i; */
/*   if (X != A && (X->n != A->n || X->m != A->m) || */
/*       X != B && (X->n != B->n || X->m != B->m))  */
/*     Error("Subtract: matrices not of same order"); */

/*   for (i = 0; i < X->n*X->m; i++) */
/*     X->a[i] = A->a[i] - B->a[i]; */
/*   return X; */
/* } */

/* Matrix *Matrix_Multiply(Matrix *X, Matrix *A, Matrix *B) */
/* { */
/*   int i, j, k; */
/*   double *a = doublevector(A->n*B->m); */

/*   if (A != B && A->m != B->n) Error("A and B cannot my multiplied"); */
/*   if (X->n != A->n && X->m != B->m) Error("X is wrong order"); */

/*   for (i = 0; i < A->n; i++) */
/*     for (j = 0; j < B->m; j++) */
/*       for (k = 0; k < A->m; k++) */
/* 	a[i*X->m+j] += A->a[i*A->m+k]*B->a[k*B->m+j]; */

/*   for (i = 0; i < X->n*X->m; i++) */
/*     X->a[i] = a[i]; */

/*   free(a); */
/*   return X; */
/* } */

/* Matrix *Matrix_Cholesky(Matrix *X, Matrix *A) */
/* { */
/*   // A must be at least lower triangular */
/*   if (X != A && (X->n != A->n || X->m != A->m))  */
/*     Error("Cholesky: matrices not of same order"); */

/*   Matrix_Copy(X, A); */
/*   if (LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', X->n, X->a, X->n)) { */
/*     fprintf(stderr, "cholesky"); */
/*     return NULL; */
/*   } */
/*   return X; */
/* } */

/* Matrix *Matrix_Invert(Matrix *X, Matrix *A) */
/* { */
/*   if (X != A && (X->n != A->n || X->m != A->m))  */
/*     Error("Invert: matrices not of same order"); */

/*   Matrix_Copy(X, A); */
/*   if (Matrix_Cholesky(X, X) != NULL && */
/*       !LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'L', X->n, X->a, X->n)) */
/*     return X; */
/*   else { */
/*     fprintf(stderr, "invert"); */
/*     return NULL; */
/*   } */
/* } */
 
/* double Matrix_Log_Determinant(Matrix *A, int n) */
/* { */
/*   int i; */
/*   double logdetA = 0; */
/*   Matrix X; */
  
/*   if (A->n != A->m) Error("det: matrix is not square"); */
/*   n = MIN(n, A->n); */
/*   Matrix_Init(&X, "det", n, n); */
/*   Matrix_Sub_Copy(&X, A, n); */

/*   if (Matrix_Cholesky(&X, &X) != NULL) */
/*     for (i = 0; i < n; i++) */
/*       logdetA += 2.*log(X.a[i*n+i]); */

/*   Matrix_Destroy(&X); */
/*   return logdetA; */
/* } */

/* void Matrix_Print(Matrix *A) */
/* { */
/*   int i, j; */
/*   printf("name=%s\n", A->name); */
/*   printf("rows=%d columns=%d\n", A->n, A->m); */
/*   for (i = 0; i < A->n; i++) { */
/*     for (j = 0; j < A->m; j++) */
/*       printf("% 7.1e ", A->a[i*A->m+j]); */
/*     printf("\n"); */
/*   } */
/* } */

double identity(double x)
{
  return x;
}

double logit(double x)
{
  return log(x/(1.-x));
}

double invlogit(double x)
{
  return 1./(1.+exp(-x));
}


static double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
/******************************************************************************/

static double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}

static double r8_error_f ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ERROR_F evaluates the error function ERF.

  Discussion:

    Since some compilers already supply a routine named ERF which evaluates
    the error function, this routine has been given a distinct, if
    somewhat unnatural, name.

    The function is defined by:

      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.

    Properties of the function include:

      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
                               ERF(0) =           0.0;
                               ERF(0.476936...) = 0.5;
      Limit ( X -> +Infinity ) ERF(X) =          +1.0.

      0.5 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2006

  Author:

    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.

  Reference:

    William Cody,
    "Rational Chebyshev Approximations for the Error Function",
    Mathematics of Computation,
    1969, pages 631-638.

  Parameters:

    Input, double X, the argument of the error function.

    Output, double R8_ERROR_F, the value of the error function.
*/
{
  double a[5] = {
    3.16112374387056560,
    1.13864154151050156E+02,
    3.77485237685302021E+02,
    3.20937758913846947E+03,
    1.85777706184603153E-01 };
  double b[4] = {
    2.36012909523441209E+01,
    2.44024637934444173E+02,
    1.28261652607737228E+03,
    2.84423683343917062E+03 };
  double c[9] = {
    5.64188496988670089E-01,
    8.88314979438837594,
    6.61191906371416295E+01,
    2.98635138197400131E+02,
    8.81952221241769090E+02,
    1.71204761263407058E+03,
    2.05107837782607147E+03,
    1.23033935479799725E+03,
    2.15311535474403846E-08 };
  double d[8] = {
    1.57449261107098347E+01,
    1.17693950891312499E+02,
    5.37181101862009858E+02,
    1.62138957456669019E+03,
    3.29079923573345963E+03,
    4.36261909014324716E+03,
    3.43936767414372164E+03,
    1.23033935480374942E+03 };
  double del;
  double erfxd;
  int i;
  double p[6] = {
    3.05326634961232344E-01,
    3.60344899949804439E-01,
    1.25781726111229246E-01,
    1.60837851487422766E-02,
    6.58749161529837803E-04,
    1.63153871373020978E-02 };
  double q[5] = {
    2.56852019228982242,
    1.87295284992346047,
    5.27905102951428412E-01,
    6.05183413124413191E-02,
    2.33520497626869185E-03 };
  double sqrpi = 0.56418958354775628695;
  double thresh = 0.46875;
  double xabs;
  double xbig = 26.543;
  double xden;
  double xnum;
  double xsmall = 1.11E-16;
  double xsq;

  xabs = fabs ( ( x ) );
/*
  Evaluate ERF(X) for |X| <= 0.46875.
*/
  if ( xabs <= thresh )
  {
    if ( xsmall < xabs )
    {
      xsq = xabs * xabs;
    }
    else
    {
      xsq = 0.0;
    }

    xnum = a[4] * xsq;
    xden = xsq;

    for ( i = 0; i < 3; i++ )
    {
      xnum = ( xnum + a[i] ) * xsq;
      xden = ( xden + b[i] ) * xsq;
    }

    erfxd = x * ( xnum + a[3] ) / ( xden + b[3] );
  }
/*
  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
*/
  else if ( xabs <= 4.0 )
  {
    xnum = c[8] * xabs;
    xden = xabs;
    for ( i = 0; i < 7; i++ )
    {
      xnum = ( xnum + c[i] ) * xabs;
      xden = ( xden + d[i] ) * xabs;
    }

    erfxd = ( xnum + c[7] ) / ( xden + d[7] );
    xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
    del = ( xabs - xsq ) * ( xabs + xsq );
    erfxd = exp ( - xsq * xsq ) * exp ( -del ) * erfxd;

    erfxd = ( 0.5 - erfxd ) + 0.5;

    if ( x < 0.0 )
    {
      erfxd = -erfxd;
    }
  }
/*
  Evaluate ERFC(X) for 4.0 < |X|.
*/
  else
  {
    if ( xbig <= xabs )
    {
      if ( 0.0 < x )
      {
        erfxd = 1.0;
      }
      else
      {
        erfxd = - 1.0;
      }
    }
    else
    {
      xsq = 1.0 / ( xabs * xabs );

      xnum = p[5] * xsq;
      xden = xsq;

      for ( i = 0; i < 4; i++ )
      {
        xnum = ( xnum + p[i] ) * xsq;
        xden = ( xden + q[i] ) * xsq;
      }

      erfxd = xsq * ( xnum + p[4] ) / ( xden + q[4] );
      erfxd = ( sqrpi - erfxd ) / xabs;
      xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
      del = ( xabs - xsq ) * ( xabs + xsq );
      erfxd = exp ( - xsq * xsq ) * exp ( - del ) * erfxd;

      erfxd = ( 0.5 - erfxd ) + 0.5;

      if ( x < 0.0 )
      {
        erfxd = -erfxd;
      }
    }
  }

  return erfxd;
}

static double r8_modp ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MODP returns the nonnegative remainder of R8 division.

  Discussion:

    If
      REM = R8_MODP ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360.0) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.

  Example:

        I         J     MOD R8_MODP  R8_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number to be divided.

    Input, double Y, the number that divides X.

    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
*/
{
  double value;

  if ( y == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_MODP - Fatal error!\n" );
    fprintf ( stderr, "  R8_MODP ( X, Y ) called with Y = %g\n", y );
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + fabs ( y );
  }

  return value;
}

double von_mises_cdf ( double x, double a, double b )

/******************************************************************************/
/*
  Purpose:

    VON_MISES_CDF evaluates the von Mises CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2006

  Author:

    Original FORTRAN77 version by Geoffrey Hill
    C version by John Burkardt

  Reference:

    Geoffrey Hill,
    ACM TOMS Algorithm 518,
    Incomplete Bessel Function I0: The von Mises Distribution,
    ACM Transactions on Mathematical Software,
    Volume 3, Number 3, September 1977, pages 279-284.

    Kanti Mardia, Peter Jupp,
    Directional Statistics,
    Wiley, 2000, QA276.M335

  Parameters:

    Input, double X, the argument of the CDF.
    A - PI <= X <= A + PI.

    Input, double A, B, the parameters of the PDF.
    -PI <= A <= PI,
    0.0 < B.

    Output, double VON_MISES_CDF, the value of the CDF.
*/
{
  double a1 = 12.0;
  double a2 = 0.8;
  double a3 = 8.0;
  double a4 = 1.0;
  double arg;
  double c;
  double c1 = 56.0;
  double cdf;
  double ck = 10.5;
  double cn;
  double erfx;
  int ip;
  int n;
  double p;
  const double r8_pi = 3.14159265358979323;
  double r;
  double s;
  double sn;
  double u;
  double v;
  double y;
  double z;
/*
  We expect -PI <= X - A <= PI.
*/
  if ( x - a <= - r8_pi )
  {
    cdf = 0.0;
    return cdf;
  }

  if ( r8_pi <= x - a )
  {
    cdf = 1.0;
    return cdf;
  }
/*
  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
*/
  z = b;

  u = r8_modp ( x - a + r8_pi, 2.0 * r8_pi );

  if ( u < 0.0 )
  {
    u = u + 2.0 * r8_pi;
  }

  y = u - r8_pi;
/*
  For small B, sum IP terms by backwards recursion.
*/
  if ( z <= ck )
  {
    v = 0.0;

    if ( 0.0 < z )
    {
      ip = ( int ) ( z * a2 - a3 / ( z + a4 ) + a1 );
      p = ( double ) ( ip );
      s = sin ( y );
      c = cos ( y );
      y = p * y;
      sn = sin ( y );
      cn = cos ( y );
      r = 0.0;
      z = 2.0 / z;

      for ( n = 2; n <= ip; n++ )
      {
        p = p - 1.0;
        y = sn;
        sn = sn * c - cn * s;
        cn = cn * c + y * s;
        r = 1.0 / ( p * z + r );
        v = ( sn / p + v ) * r;
      }
    }
    cdf = ( u * 0.5 + v ) / r8_pi;
  }
/*
  For large B, compute the normal approximation and left tail.
*/
  else
  {
    c = 24.0 * z;
    v = c - c1;
    r = sqrt ( ( 54.0 / ( 347.0 / v + 26.0 - c ) - 6.0 + c ) / 12.0 );
    z = sin ( 0.5 * y ) * r;
    s = 2.0 * z * z;
    v = v - s + 3.0;
    y = ( c - s - s - 16.0 ) / 3.0;
    y = ( ( s + 1.75 ) * s + 83.5 ) / v - y;
    arg = z * ( 1.0 - s / y / y );
    erfx = r8_error_f ( arg );
    cdf = 0.5 * erfx + 0.5;
  }
  cdf = r8_max ( cdf, 0.0 );
  cdf = r8_min ( cdf, 1.0 );

  return cdf;
}
/******************************************************************************/

double von_mises_cdf_inv ( double cdf, double a, double b )

/******************************************************************************/
/*
  Purpose:

    VON_MISES_CDF_INV inverts the von Mises CDF.

  Discussion:

    A simple bisection method is used on the interval [ A - PI, A + PI ].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double CDF, the value of the CDF.

    Input, double A, B, the parameters of the PDF.
    -PI <= A <= PI,
    0.0 < B.

    Output, double VON_MISES_CDF_INV, the corresponding argument of the CDF.
    A - PI <= X <= A + PI.
*/
{
  double cdf1;
  double cdf3;
  int it;
  int it_max = 100;
  const double r8_pi = 3.14159265358979323;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "VON_MISES_CDF_INV - Fatal error!\n" );
    fprintf ( stderr, "  CDF < 0 or 1 < CDF.\n" );
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a - r8_pi;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = a + r8_pi;
    return x;
  }
  x1 = a - r8_pi;
  cdf1 = 0.0;

  x2 = a + r8_pi;
/*
  Now use bisection.
*/
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = von_mises_cdf ( x3, a, b );

    if ( fabs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "VON_MISES_CDF_INV - Fatal error!\n" );
      fprintf ( stderr, "  Iteration limit exceeded.\n" );
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
    }
  }

  return x;
}
/******************************************************************************/
