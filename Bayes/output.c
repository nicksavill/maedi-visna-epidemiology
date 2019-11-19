#include <bayes.h>

extern int nt1, nt2, nr1, nr2;

int Samples(struct simulation *Sim, struct data *Data, double *****samples);
void ReadBestParameters(Simulation *Sim, struct parameterset *Q,
            int nt, int nr, double *p, char* infile);
int HyperSamples(struct simulation *Sim, struct data *Data, double ****samples);
void RunModel(struct data *Data, struct variables *V, int i, int maxs, int nss, 
          int nt, int nr, double *theta, 
          struct data_individual *Ind, 
          double ***residual, double ***preddist, double ***simdist,
          double ***fitdist, double ***lnL, gsl_rng *stream);

int Samples(struct simulation *Sim, struct data *Data, double *****samples)
{
  int i, k, ns, nt, nr, npar, *use;
  double *p;
  char filename[1000];
  struct parameterset *Q;
  FILE *fp;

  if (!(fp = fopen(Sim->logfile, "r")))
    Error("no posterior file");
  fclose(fp);
  if (Sim->mcmc)
    sprintf(filename,"Results/%s.post", Sim->outexpno);
  else
    sprintf(filename,"Results/%s.post", Sim->expno);
  if (!(fp = fopen(filename, "r")))
    Error("no posterior file");

  fscanf(fp, "%*s%d", &ns);
  *samples = (double****) malloc(Data->ntrt*sizeof(double***));
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    (*samples)[nt] = (double***) malloc(Data->Treatment[nt].nind*sizeof(double**));
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
      (*samples)[nt][nr] = 
    doublematrix(ns, Sim->Treatment[nt].Individual[nr].Q.parameters+3);
  }
  
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;

      p = doublevector(Q->parameters);
      Parameters(nt, nr, p, Q);
      for (k = 0; k < Q->parameters; k++) 
        if (Use(INDIVIDUAL_CONSTANT, Q->P[k].type)) 
          Q->P[k].location = p[k];
      free(p);
      
      use = integervector(Q->parameters);
      fscanf(fp, "%*s%*d%*d%d", &npar);
      for (i = 0; i < npar; i++) {
        fscanf(fp, "%d", &k);
        use[k] = TRUE;
      }

      for (i = 0; i < ns; i++) {
    for (k = 0; k < Q->parameters; k++)
      if (use[k]) {
        fscanf(fp, "%lf", &((*samples)[nt][nr][i][k]));
      }
      else if (Use(CONSTANT, Q->P[k].type) ||
           Use(INDIVIDUAL_CONSTANT, Q->P[k].type))
        (*samples)[nt][nr][i][k] = Q->P[k].location;
    for (k = Q->parameters; k < Q->parameters+3; k++)
      fscanf(fp, "%lf", &((*samples)[nt][nr][i][k]));
      }
      free(use);
    }
  fclose(fp);

  return ns;
}

void ReadBestParameters(Simulation *Sim, struct parameterset *Q,
            int nt, int nr, double *p, char* filename)
{
  int i, j, k;
  double v;
  FILE *fp;
  struct parameter *P = Q->P;
  
  for (k = 0; k < Q->parameters; k++)
    if (Use(CONSTANT, P[k].type))
      p[k] = P[k].location;

  fp = fopen(filename, "r");
  while (fscanf(fp, "%d%d%d%lf", &i, &j, &k, &v) != EOF) {
    if (!Use(CONSTANT, P[k].type) && P[k].type != NOTUSED && i == nt && j == nr)
      p[k] = v;
  }
  
  fclose(fp);
}

int HyperSamples(struct simulation *Sim, struct data *Data, double ****samples)
{
  int i, k, ns, nt, npar;
  char filename[1000];
  FILE *fp;

  if (!(fp = fopen(Sim->logfile, "r")))
    Error("no posterior file");
  fclose(fp);
  if (Sim->mcmc)
    sprintf(filename,"Results/%s.hypr", Sim->outexpno);
  else
    sprintf(filename,"Results/%s.hypr", Sim->expno);
  if (!(fp = fopen(filename, "r")))
    Error("no hyper posterior file");

  fscanf(fp, "%*s%d", &ns);
  *samples = (double***) malloc(Data->ntrt*sizeof(double**));
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    (*samples)[nt] = doublematrix(ns, Sim->HP.nhyper_parameters);
  
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    fscanf(fp, "%*s%*d%d", &npar);
    for (i = 0; i < npar; i++)
      fscanf(fp, "%*d");

    for (i = 0; i < ns; i++)
      for (k = 0; k < npar; k++)
      fscanf(fp, "%lf", &((*samples)[nt][i][k]));
    }
  fclose(fp);
  return ns;
}

void RecalculateDerivedParameters(struct simulation *Sim, struct data *Data)
{
  int i, nt, nr, k, ns;
  double ****samples, *theta, **y;
  char filename[1000];
  struct parameterset *Q;
  struct data_individual *Ind;
  FILE *fp;
  TimePoints TP;

  if (!Sim->rdp)
    return;
  else if (Sim->verbose)
    printf("Recalculating derived parameters\n");

  ns = Samples(Sim, Data, &samples);
  sprintf(filename,"Results/%s.post", Sim->outexpno);
  fp = fopen(filename, "w");

  fprintf(fp, "samples %d\n", ns);
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Data->Treatment[nt].Individual[nr];
      Q = &Sim->Treatment[nt].Individual[nr].Q;

      fprintf(fp, "individual %d %d %d ", nt, nr, Q->npar);
      for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) 
      fprintf(fp, "%d ", k);
    fprintf(fp, "\n");

      for (i = 0; i < ns; i++) {
    theta = samples[nt][nr][i];
    for (k = 0; k < Q->parameters; k++) 
      if (Use(CONSTANT, Q->P[k].type) || 
          Use(INDIVIDUAL_CONSTANT, Q->P[k].type))
        theta[k] = Q->P[k].location;

    
    TP.mode = OUTPUT;
    TP.n = Ind->np;
    TP.i = 0;
    TP.t = Ind->Y;
    // TP.t = doublevector(TP.n);
    // for (k = 0; k < TP.n; k++)
    //   TP.t[k] = Ind->Y[k][0];

    k = nsteps(Sim, Ind);
    y = doublematrix(k, 1+Sim->V.variables);
    Clear(y, k, Sim->V.variables);
    DerivedParameters(nt, nr, k, theta, y, &TP);
    // free(TP.t); 
    free(y[0]); free(y);
    
    for (k = 0; k < Q->parameters; k++) 
      if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) {
        if (Q->P[k].number == REAL)
          fprintf(fp, "%.3e ", theta[k]);
        else if (Q->P[k].number == INTEGER)
          fprintf(fp, "%ld ", lrint(theta[k]));
      }
    fprintf(fp, "%.3e %.3e %ld\n", theta[k], theta[k+1], lrint(theta[k+2]));
      }
    }
  fclose(fp);
}

void PostAnalysis(struct simulation *Sim, struct data *Data)
{
  double ****samples, **d1;
  int f, ns, nt, nr, i, k, np1 = 0, nl50, nu50, nl95, nu95, me;
  char filename[100];
  FILE *fpo;
  TimePoints TP;
  struct data_individual *Ind;
  
  if (Sim->verbose)
    printf("Post analysis\n");

  // open file for writing
  sprintf(filename,"Results/%s.func", Sim->outexpno);
  fpo = fopen(filename, "w");

  // get the posterior samples and the number of those samples
  ns = Samples(Sim, Data, &samples);

  // define quantiles
  nl50 = lrint(0.25*ns);
  nu50 = lrint(0.75*ns);
  nl95 = lrint(0.05*ns);
  nu95 = lrint(0.95*ns);
  me = lrint(0.5*ns);

  // get the list of functions to calculate
  functions_t F;
  function_list(&F);

  // run each function on the samples
  for (f = 0; f < F.n_func; f++) {
    // get the number of points for this function
    np1 = F.function_list[f](nt1, nr1, NULL, NULL, NULL, &TP);
    if (np1 < 1 || np1 > 1e8)
      continue;

    // get x for output
    double *x = doublevector(np1);
    double *y = doublevector(np1);
    F.function_list[f](nt1, nr1, x, y, samples[nt1][nr1][0], &TP);

    // matrix of stored values 
    d1 = doublematrix(np1, ns);

    //output information about this function
    fprintf(fpo, "function %d\n", f);
    fprintf(fpo, "%d\n", np1);
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
        fprintf(fpo, "%d %d\n", nt, nr);
    fprintf(fpo, "\n");

    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
        Ind = &Data->Treatment[nt].Individual[nr];

        TP.mode = OUTPUT;
        TP.n = Ind->np;
        TP.i = 0;
        TP.t = Ind->Y;
        //   TP.t = doublevector(TP.n);
        //   for (k = 0; k < TP.n; k++)
        // TP.t[k] = Ind->Y[k][0];

        // count = 0;
#       pragma omp parallel private(i, k) num_threads(Sim->threads)
        {
          double *x1 = doublevector(np1);
          double *y1 = doublevector(np1);
          int iam = omp_get_thread_num();
          int nth = omp_get_num_threads();
          int ip = ns/nth;
          int is = iam*ip;
          if (iam == nth-1) ip = ns-is;

          for (i = is; i < is+ip; i++) {
            F.function_list[f](nt, nr, x1, y1, samples[nt][nr][i], &TP);
#           pragma omp critical (covariance)
            for (k = 0; k < np1; k++)
              d1[k][i] = y1[k];
          }
          free(x1);
          free(y1);
        }
        for (k = 0; k < np1; k++)
          HeapSort(d1[k]-1, ns);
        for (k = 0; k < np1; k++)
          fprintf(fpo, "%e %.3e %.3e %.3e %.3e %.3e\n", x[k], d1[k][nl95], 
          d1[k][nu95], d1[k][nl50], d1[k][nu50], d1[k][me]);
        fprintf(fpo, "\n");
      }
      free(x);
      free(y);
      free(d1[0]);
      free(d1);
  }
}

/* static void hpd(double *x, double alpha, double *low, double *high, int n) */
/* { */
/*   int start = 0; */
/*   int end = lrint(n*(1.-alpha)); */
/*   double width, min_width = INFINITY; */

/*   while(end < n) { */
/*     width = x[end] - x[start]; */
/*     if (width < min_width) { */
/*       min_width = width; */
/*       *low = x[start]; */
/*       *high = x[end]; */
/*     } */
/*     start++; */
/*     end++; */
/*   } */
/* } */

void Marginals(struct simulation *Sim, struct data *Data)
{
  int i, k, n, nt, nr, ns;
  //  int hns, g, np = Sim->hypers_used+Sim->derived;
  double **y;
  double ****samples, ****dist, ***mle, ***imean, ***ivar;
  //  double ***hsamples, ***hdist, **himean, hivar;
  char filename[100];
  struct parameterset *Q;
  struct data_individual *Ind;
  FILE *fps;
  TimePoints TP;
  
  if (!Sim->marg) return;
  if (Sim->verbose) printf("Marginal distributions\n");
  sprintf(filename,"Results/%s.marg", Sim->outexpno);
  fps = fopen(filename, "w");

  /* hyper parameters */
  /* if (Sim->hypers_used) { */
  /*   hns = HyperSamples(Sim, Data, &hsamples); */
  /*   hdist  = (double***) malloc(Data->ntrt*sizeof(double**)); */
  /*   himean = (double**)  malloc(Data->ntrt*sizeof(double*)); */
  /*   hivar  = (double**)  malloc(Data->ntrt*sizeof(double*)); */
  /*   for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) { */
  /*     hdist[nt]  = doublematrix(np, hns); */
  /*     himean[nt] = doublevector(np); */
  /*     hivar[nt]  = doublevector(np); */
  /*   } */

  /*   for (i = 0; i < hns; i++) */
  /*     for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) */
  /*     for (k = 0; k < np; k++)  */
  /*       hdist[nt][k][i] = hsamples[nt][i][k]; */
  
  /*   for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) */
  /*     for (k = 0; k < np; k++) { */
  /*     HeapSort(hdist[nt][k]-1, hns); */
  /*     meanvar(hdist[nt][k], hns, &himean[nt][k], &hivar[nt][k]); */
  /*     } */

  /*   /\* output hyper parameters *\/ */
  /*   fprintf(fps, "HYPERS\ntreatment parameter names\n"); */
  /*   for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) { */
  /*     HS = &Sim->H; */
  /*     for (i = 0; i < HS->npar; i++) { */
  /*     k = HS->par_indices[i]; */
  /*     fprintf(fps, "%d %d \"%s\"\n", nt, k, HS->name[i]); */
  /*     } */
  /*     HT = &Sim->Treatment[nt].H; */
  /*     for (i = 0; i < HT->npar; i++) { */
  /*     k = HT->par_indices[i]; */
  /*     fprintf(fps, "%d %d \"%s\"\n", nt, k, HT->name[i]); */
  /*     } */
  /*   } */

  /*   fprintf(fps, "treatment parameters median (50%% CI) (95%% CI) mean sd\n"); */
  /*   for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) { */
  /*     HS = &Sim->H; */
  /*     for (i = 0; i < HS->npar; i++) { */
  /*     k = HS->par_indices[i]; */
  /*     fprintf(fps, "%d %d %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",  */
  /*         nt, k, */
  /*         hdist[nt][i][lrint(0.5*hns)], */
  /*         hdist[nt][i][lrint(0.25*hns)], */
  /*         hdist[nt][i][lrint(0.75*hns)], */
  /*         hdist[nt][i][lrint(0.025*hns)], */
  /*         hdist[nt][i][lrint(0.975*hns)], */
  /*         himean[nt][i], sqrt(hivar[nt][i])); */
  /*     } */
  /*     HT = &Sim->Treatment[nt].H; */
  /*     for (i = 0; i < HT->npar; i++) { */
  /*     k = HT->par_indices[i]; */
  /*     fprintf(fps, "%d %d %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",  */
  /*         nt, k, */
  /*         hdist[nt][i][lrint(0.5*hns)], */
  /*         hdist[nt][i][lrint(0.25*hns)], */
  /*         hdist[nt][i][lrint(0.75*hns)], */
  /*         hdist[nt][i][lrint(0.025*hns)], */
  /*         hdist[nt][i][lrint(0.975*hns)], */
  /*         himean[nt][i], sqrt(hivar[nt][i])); */
  /*     } */
  /*   } */
  /* } */

  /* individual level parameters */
  ns  = Samples(Sim, Data, &samples);
  dist  = (double****) malloc(Data->ntrt*sizeof(double***));
  mle   = (double***)  malloc(Data->ntrt*sizeof(double**));
  imean = (double***)  malloc(Data->ntrt*sizeof(double**));
  ivar  = (double***)  malloc(Data->ntrt*sizeof(double**));
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    dist[nt]  = (double***) malloc(Data->Treatment[nt].nind*sizeof(double**));
    mle[nt]   = (double**)  malloc(Data->Treatment[nt].nind*sizeof(double*));
    imean[nt] = (double**)  malloc(Data->Treatment[nt].nind*sizeof(double*));
    ivar[nt]  = (double**)  malloc(Data->Treatment[nt].nind*sizeof(double*));

    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      dist[nt][nr]  = doublematrix(Q->parameters, ns);
      mle[nt][nr]   = doublevector(Q->parameters);
      imean[nt][nr] = doublevector(Q->parameters);
      ivar[nt][nr]  = doublevector(Q->parameters);
    }
  }

  for (i = 0; i < ns; i++)
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
    Q = &Sim->Treatment[nt].Individual[nr].Q;
        for (k = 0; k < Q->parameters; k++) 
      if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type))
        dist[nt][nr][k][i] = samples[nt][nr][i][k];
      }
  
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) {
      HeapSort(dist[nt][nr][k]-1, ns);
      if (Q->P[k].circular)
        circmeanvar(dist[nt][nr][k], ns, &imean[nt][nr][k], 
            &ivar[nt][nr][k]);
      else
        meanvar(dist[nt][nr][k], ns, &imean[nt][nr][k], 
            &ivar[nt][nr][k]);
    }
    }

  /* read in MLE parameters and calculate derived parameters */
  sprintf(filename, "Results/%s.best", Sim->outexpno);
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      Ind = &Data->Treatment[nt].Individual[nr];
      ReadBestParameters(Sim, Q, nt, nr, mle[nt][nr], filename);

      TP.mode = OUTPUT;
      TP.n = Ind->np;
      TP.i = 0;
      TP.t = Ind->Y;
    //   TP.t = doublevector(TP.n);
    //   for (k = 0; k < TP.n; k++)
    // TP.t[k] = Ind->Y[k][0];

      n = nsteps(Sim, Ind);
      y = doublematrix(n, 1+Sim->V.variables);
      Clear(y, n, Sim->V.variables);
      DerivedParameters(nt, nr, n, mle[nt][nr], y, &TP);
      // free(TP.t); 
      free(y[0]); free(y);
    }

  /* output individual level parameters */
  fprintf(fps, "INDIVIDUALS\ntreatment individual parameter name\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type))
      fprintf(fps, "%d %d %d \"%s\"\n", nt, nr, k, Q->P[k].name);
    }

  fprintf(fps, "parameter treatment individual median (50%% CI) (95%% CI) mean sd "
      "MAP\n");
  for (k = 0; k < Sim->maxnp; k++)
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
    Q = &Sim->Treatment[nt].Individual[nr].Q;
    if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) {
      fprintf(fps, "%d %d %d %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", 
          k, nt, nr,
          dist[nt][nr][k][lrint(0.5*ns)],
          dist[nt][nr][k][lrint(0.25*ns)],
          dist[nt][nr][k][lrint(0.75*ns)],
          dist[nt][nr][k][lrint(0.025*ns)],
          dist[nt][nr][k][lrint(0.975*ns)],
          // dist[nt][nr][k][lrint(0.005*ns)],
          // dist[nt][nr][k][lrint(0.995*ns)],
          imean[nt][nr][k], sqrt(ivar[nt][nr][k]), mle[nt][nr][k]);
    }
      }
  fclose(fps);
}

void RunModel(struct data *Data, struct variables *V, int i, int maxs, int nss, 
          int nt, int nr, double *theta, struct data_individual *Ind, 
          double ***residual, double ***lnL, double ***preddist, double ***simdist,
          double ***fitdist, gsl_rng *stream)
{
  int k, l, status;
  int *t_index = integervector(Ind->np); 
  uint *truevector = uintegervector(Data->nvar);
  double **y = doublematrix(nss, V->variables);
  double *Y  = doublevector(max(Data->nvar, V->variables));
  double **O2 = doublematrix(Ind->np, V->variables);
  double *O  = doublevector(V->variables);
  double *Oo = doublevector(V->variables);
  TimePoints TP;

  TP.mode = OUTPUT;
  TP.n = Ind->np;
  TP.i = 0;
  TP.t = Ind->Y;
  // TP.t = doublevector(TP.n);
  // for (k = 0; k < TP.n; k++)
  //   TP.t[k] = Ind->Y[k][0];
  
  for (l = 0; l < Data->nvar; l++)
    truevector[l] = Measured;

  Clear(y, nss, V->variables);
  if ((status = Model(nt, nr, nss, theta, y, &TP))) {
    fprintf(stderr, "Model solution has changed. Rerun full MCMC simulation: %d %d\n",
        nt, nr);
    fprintf(stderr, "Error # %d\n", status);
    Error("");
  }
  Time_Index(Ind, nss, y, t_index);

  for (k = 0; k < Ind->np; k++) {
    SetNAN(O, V->variables);
    Residual(nt, nr, O, y[t_index[k]], theta, Ind->Y[k], Ind->value[k]);
    for (l = 1; l < V->variables; l++)
      residual[l][k][i] = O[l];
  }

  for (k = 0; k < Ind->np; k++) {
    SetNAN(O, V->variables);
    SetNAN(Y, V->variables);
    SimulateData(nt, nr, Y, y[t_index[k]], theta, Ind->Y[k], Ind->value[k], stream);
    OutputData(nt, nr, O, Oo, Y, truevector, k);
    for (l = 1; l < V->variables; l++)
      simdist[l][k][i] = O[l];
  }

  // for (k = 0; k < Ind->np; k++) {
  //   SetNAN(O, V->variables);
  //   WAIC(nt, nr, O, y[t_index[k]], theta, Ind->Y[k], Ind->value[k]);
  //   for (l = 1; l < V->variables; l++)
  //     lnL[l][k][i] = O[l];
  // }
  for (k = 0; k < Ind->np; k++)
    for (l = 1; l < V->variables; l++)
      SetNAN(O2[k], V->variables);
  WAIC(nt, nr, O2, y, theta, &TP);
  for (k = 0; k < Ind->np; k++)
    for (l = 1; l < V->variables; l++)
      lnL[l][k][i] = O2[k][l];
  

  if (i < maxs) {
    for (k = 0; k < nss; k++) {
      SetNAN(O, V->variables);
      OutputModel(nt, nr, O, y[k]);
      for (l = 1; l < V->variables; l++)
    fitdist[l][k][i] = O[l];
    }
    
    for (k = 0; k < nss; k++) {
      SetNAN(O, V->variables);
      PredictData(nt, nr, Y, y[k], theta, stream);
      OutputData(nt, nr, O, Oo, Y, truevector, k);
      for (l = 1; l < V->variables; l++)
    preddist[l][k][i] = O[l];
    }
    
  }
  free(t_index); free(O); free(Y); free(y[0]); free(y); free(truevector);
  free(Oo); 
  free(O2[0]); free(O2);
  // free(TP.t);
}

void PosteriorPrediction(struct simulation *Sim, struct data *Data)
{
  int i, k, nt, nr, ns, l, nss, *t_index, count;
  int nl95, nu95, nl50, nu50, median, maxs, firsttime;
  double R, sumR, sumR2, *Ysim, *****Yind, ***residual;
  double *****lnL, meansim, varsim;
  double **theta, ****samples, **y, ***fitdist, ***preddist, ***simdist;
  double *O = doublevector(Sim->V.variables);
  double *Oo = doublevector(Sim->V.variables);
  char filename[100];
  struct data_individual *Ind;
  struct data_treatment *T;
  struct variables *V = &Sim->V;
  struct parameterset *Q;
  //  gsl_rng *stream[Sim->threads];
  FILE *fps, *fpw;
  TimePoints TP;
  long double meanlnL, varlnL;//, *tmplnL;
  // int countlnL;

  if (!Sim->popd) return;
  ns = Samples(Sim, Data, &samples);
  if (Sim->verbose) printf("Posterior prediction\n");
  
  Yind = (double*****) malloc(Data->ntrt*sizeof(double****));
  lnL = (double*****) malloc(Data->ntrt*sizeof(double****));
  // tmplnL = (long double*) malloc(ns*sizeof(long double));

  sprintf(filename,"Results/%s.stat", Sim->outexpno);
  fps = fopen(filename, "w");
  sprintf(filename,"Results/%s.waic", Sim->outexpno);
  fpw = fopen(filename, "w");
  sprintf(filename, "Results/%s.best", Sim->outexpno);

  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    Yind[nt] = (double****) malloc(Data->Treatment[nt].nind*sizeof(double***));
    lnL[nt] = (double****) malloc(Data->Treatment[nt].nind*sizeof(double***));
    T = &Data->Treatment[nt];
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      Yind[nt][nr] = doublematrix3(V->variables, Ind->np, ns);
      lnL[nt][nr] = doublematrix3(V->variables, Ind->np, ns);

      nss = nsteps(Sim, Ind);
      maxs = MIN(ns, lrint(pow(2, 26)/(V->variables*nss)));
      t_index = integervector(Ind->np); 
      y = doublematrix(nss, V->variables);
      preddist = doublematrix3(V->variables, nss, maxs);
      simdist = doublematrix3(V->variables, Ind->np, maxs);
      fitdist = doublematrix3(V->variables, nss, maxs);
      Ysim = doublevector(V->variables);
      theta = samples[nt][nr];
      residual = Yind[nt][nr];

      nl50 = lrint(0.25*maxs);
      nu50 = lrint(0.75*maxs);
      nl95 = lrint(0.025*maxs);
      nu95 = lrint(0.975*maxs);
      median = lrint(0.5*maxs);
      
#     pragma omp parallel private(i) num_threads(Sim->threads)
      {
        int iam = omp_get_thread_num();
        int nth = omp_get_num_threads();
        int ip = ns/nth;
        int is = iam*ip;
        if (iam == nth-1) ip = ns-is;

        gsl_rng *stream = gsl_rng_alloc(gsl_rng_taus);    
        gsl_rng_set(stream, Sim->seed+iam);

        for (i = is; i < is+ip; i++)
          RunModel(Data, V, i, maxs, nss, nt, nr, theta[i],
               Ind, residual, lnL[nt][nr], 
               preddist, simdist, fitdist, 
               stream);
        gsl_rng_free(stream);
      }
      TP.mode = OUTPUT;
      TP.n = Ind->np;
      TP.i = 0;
      TP.t = Ind->Y;
      // TP.t = doublevector(TP.n);
      // for (k = 0; k < TP.n; k++)
      //   TP.t[k] = Ind->Y[k][0];
      
      ReadBestParameters(Sim, Q, nt, nr, theta[0], filename);
      Clear(y, nss, V->variables);
      Model(nt, nr, nss, theta[0], y, &TP);
      // Time_Index(Ind, nss, y, t_index);
      // free(TP.t);
      
      fprintf(fps, "solutions %d %d\n", nt, nr);
      for (l = 1; l < V->variables; l++) {
        fprintf(fps, "%d\n", l);
        for (k = 0; k < nss; k++) {
          SetNAN(O, V->variables);
          OutputModel(nt, nr, O, y[k]);
          HeapSort(fitdist[l][k]-1, maxs);
          meanvar(fitdist[l][k], maxs, &meansim, &varsim);
          fprintf(fps,"%.3e %.3e %.3e %.3e %.3e %.3e %.3e ", 
              y[k][0]*V->scaletime, O[l], fitdist[l][k][median],
              // y[k][0]*V->scaletime, O[l], meansim,
              fitdist[l][k][nl50], fitdist[l][k][nu50],
              fitdist[l][k][nl95], fitdist[l][k][nu95]);
          HeapSort(preddist[l][k]-1, maxs);
          if (!isnan(preddist[l][k][maxs-1]))
            fprintf(fps,"%.3e %.3e %.3e %.3e",
                preddist[l][k][nl50], preddist[l][k][nu50],
                preddist[l][k][nl95], preddist[l][k][nu95]);
          fprintf(fps, "\n");
        }
      }
      fflush(fps);

      fprintf(fps, "data %d %d\n", nt, nr);
      for (l = 1; l < V->variables; l++) {
        firsttime = TRUE;
        for (k = 0; k < Ind->np; k++) {
          SetNAN(O, V->variables);
          SetNAN(Oo, V->variables);
          OutputData(nt, nr, O, Oo, Ind->Y[k], Ind->value[k], k);
          if (!isnan(O[l])) {
            if (firsttime) {
              fprintf(fps, "%d\n", l);
              firsttime = FALSE;
            }
            fprintf(fps, "%.3e %.3e %.3e\n", Ind->Y[k][0]*V->scaletime, O[l],
                Oo[l]);
          }
        }
      }
      fflush(fps);

      fprintf(fps, "standardised residuals within individual %d %d\n", nt, nr);
      for (l = 1; l < V->variables; l++) {
        firsttime = TRUE;
        for (k = 0; k < Ind->np; k++) {

          // residual based on posterior predictive distribution
          // observed data
          // SetNAN(O, V->variables);
          // OutputData(nt, nr, O, Oo, Ind->Y[k], Ind->value[k], k);
          // // expectation and variance
          // HeapSort(simdist[l][k]-1, maxs);
          // meanvar(simdist[l][k], maxs, &meansim, &varsim);
          // if (varsim == 0) {
          //   // printf("%d %d %d %d %f varsim=%f (O=%f E=%f)\n", nt, nr, l, k,
          //   //    Ind->Y[k][l], varsim, O[l], meansim);
          //   R = 0;
          // }
          // else
          //   /* R = (O[l]-meansim)/sqrt(varsim); */
          //   R = (O[l]-simdist[l][k][median])/sqrt(varsim);

          // sumR and meanlnL are not needed 
          // residual based on 95%CI of posterior distribution
          sumR = 0;
          sumR2 = 0;
          count = 0;
          // countlnL = 0;
          for (i = 0; i < ns; i++) {

            if (!isnan(residual[l][k][i])) {
              sumR += residual[l][k][i];
              sumR2 += residual[l][k][i]*residual[l][k][i];
              count++;
            }

            // if (!isnan(lnL[nt][nr][l][k][i])) {
            //   tmplnL[countlnL] = lnL[nt][nr][l][k][i];
            //   countlnL++;
            // }
          }
      
      // if (countlnL)
      //   meanvarl(tmplnL, countlnL, &meanlnL, &varlnL);
        
      if (firsttime) {
        fprintf(fps, "%d\n", l);
        firsttime = FALSE;
      }
      sumR /= (double)count;
      sumR2 /= (double)count;
      fprintf(fps, "%.3e %.3e %.3e %.3e %.3Le\n",
          y[t_index[k]][0]*V->scaletime, R, sumR, sumR2, meanlnL);
    }
      }

      fflush(fps); fflush(fpw);
      free(y[0]); free(y); free(t_index); free(Ysim);
      free(fitdist[0][0]); free(fitdist[0]); free(fitdist); 
      free(preddist[0][0]); free(preddist[0]); free(preddist); 
      free(simdist[0][0]); free(simdist[0]); free(simdist); 
      }
  }

  long double lppd, pwaic1, pwaic2, meanL, varL, *tmp1,
    *tmp2;         
  long double Vlppd, Vpwaic1;
  tmp1 = (long double*) malloc(ns*sizeof(long double));
  tmp2 = (long double*) malloc(ns*sizeof(long double));
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    T = &Data->Treatment[nt];
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      lppd = pwaic1 = pwaic2 = 0;
      Vlppd = Vpwaic1 = 0;
      for (l = 1; l < V->variables; l++) {
    for (k = 0; k < Ind->np; k++) {
      count = 0;
      for (i = 0; i < ns; i++)
        if (!isnan(lnL[nt][nr][l][k][i])) {
          tmp1[count] = expl(lnL[nt][nr][l][k][i]);
          tmp2[count] = lnL[nt][nr][l][k][i];
          count++;
        }
      if (count) {
        wmeanvar(tmp1, count, &meanL,   &varL);
        meanvarl(tmp2, count, &meanlnL, &varlnL);
        //        varlnL *= meanlnL*meanlnL;
        //        varL is approx. meanL^2 * varL / meanL^2 = varL
        varL /= (double)count;
        varlnL /= (double)count;
        lppd += logl(meanL);
        Vlppd += varL;
        /* if (isnan(varL)) { */
        /*   /\* for (i = 0; i < count; i++) *\/ */
        /*   /\*     printf("%d %Le %Le\n", i, tmp2[i], tmp1[i]); *\/ */
        /*   printf("varL: %d %d %d, %Le\n", nt, nr, count, meanL); */
        /*   printf("%Le %Le %Le %Le\n", varL, meanlnL, varlnL, lppd); */
        /*   exit(0); */
        /* } */
        /* if (isnan(varlnL)) printf("varlnL: %d %d %d, %Lf\n", nt, nr, */
        /*                   count, meanlnL); */
        pwaic1 += 2.*(logl(meanL)-meanlnL);
        Vpwaic1 += 4.*(varL+varlnL);
        pwaic2 += varlnL;
      }
    }
      }
      fprintf(fpw, "%d %d %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le\n", nt, nr, lppd,
          pwaic1, pwaic2, -2.*(lppd-pwaic1), -2.*(lppd-pwaic2),
          4.*(Vlppd+Vpwaic1));
    }
  }

  free(O);
  fclose(fps);
  fclose(fpw);
}

void SimulatedData(struct simulation *Sim, struct data *Data)
{
  // uses pars_file to input parameter values to simulate data from
  int nt, nr, l, k, nss, *t_index;
  double **y, *Ysim, *theta;
  char filename[100];
  FILE *fpsim;
  struct data_individual *Ind;
  struct data_treatment *T;
  struct parameterset *Q;
  gsl_rng *stream;
  TimePoints TP;
  
  if (!Sim->sim) return;
  if (Sim->verbose) printf("Simulated data\n");

  sprintf(filename,"Results/%s.sim", Sim->outexpno);
  fpsim = fopen(filename, "w");

  stream = gsl_rng_alloc(gsl_rng_taus);    
  gsl_rng_set(stream, Sim->seed);

  fprintf(fpsim, "Variables %d\n", Data->nvar-1);

  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    T = &Data->Treatment[nt];
    fprintf(fpsim, "Treatment ");
    for (l = 1; l < Data->nvar; l++)
      if (T->usevar[l])
        fprintf(fpsim, "1 ");
      else
        fprintf(fpsim, "0 ");
    fprintf(fpsim, "\n");
    
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      fprintf(fpsim, "Individual\n");
      Ind = &T->Individual[nr];
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      t_index = integervector(Ind->np);
      nss = nsteps(Sim, Ind);
      y = doublematrix(nss, 1+Sim->V.variables);
      Ysim = doublevector(1+Sim->V.variables);
      theta = doublevector(Sim->Treatment[nt].Individual[nr].Q.parameters);

      TP.mode = OUTPUT;
      TP.n = Ind->np;
      TP.i = 0;
      TP.t = Ind->Y;
    //   TP.t = doublevector(TP.n);
    //   for (j = 0; j < TP.n; j++)
    // TP.t[j] = Ind->Y[j][0];
      
      ReadBestParameters(Sim, Q, nt, nr, theta, Sim->parsfile);
      Clear(y, nss, Sim->V.variables);
      Model(nt, nr, nss, theta, y, &TP);
      Time_Index(Ind, nss, y, t_index);

      for (k = 0; k < Ind->np; k++) {
    SimulateData(nt, nr, Ysim, y[t_index[k]], theta, Ind->Y[k], 
             Ind->value[k], stream);
        fprintf(fpsim, "\t%.3e", Ind->Y[k][0]);
        for (l = 1; l < Data->nvar; l++) if (T->usevar[l]) {
      if (Ind->value[k][l] == Below_Detection)
        fprintf(fpsim, "\tBD");
      else if (Ind->value[k][l] == Above_Detection)
        fprintf(fpsim, "\tAD");
      else if (Ind->value[k][l] == Measured)
        fprintf(fpsim, "\t%.3e", Ysim[l]);
      else
        fprintf(fpsim, "\tNA");
    }
        fprintf(fpsim, "\n");
      }
      
      free(theta); free(y[0]); free(y); free(t_index); free(Ysim); 
      // free(TP.t);
    }
  }
  gsl_rng_free(stream);
}
