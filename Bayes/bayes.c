#include <bayes.h>

int nt1 = 0, nt2 = -1, nr1 = 0, nr2 = -1;
char parameterfilename[50] = "pars";
char parameterlogfilename[50] = "/dev/null";

int main(int argc, char *argv[])
{
  int nt, nr;
  time_t t;
  char treatments[1000] = "\0", individuals[1000] = "\0", ignoreinds[1000] = "\0";
  char timestart[100] = "NA", timeend[100] = "NA";
  char backup[1000];
  struct data Data;  
  Simulation Sim;
  FILE *fp;

  /* -------- ignore hangup signals (eg closing terminals and internet connections */
  signal(SIGHUP, SIG_IGN); 

  /* ------------------------------------------ initialise default McMC parameters */
  Default(&Sim);

  /* ----------------------------------------------- initialise parameter log file */

  if (argc > 1) strcpy(parameterfilename, argv[1]);
  ReadParam("int", &Sim.verbose, "verbose");
  ReadParam("int", &Sim.mcmc, "mcmc");
  ReadParam("int", &Sim.rdp,  "derived");
  ReadParam("char*", Sim.expno, "expno");
  strcpy(Sim.outexpno, Sim.expno);
  ReadParam("char*", Sim.outexpno, "outexpno");
  sprintf(Sim.logfile,   "Results/%s.log",  Sim.outexpno);

  sprintf(backup, "cp %s %s.bak", Sim.logfile, Sim.logfile);
  system(backup);

  if (Sim.mcmc || Sim.rdp) {
    sprintf(parameterlogfilename, "Results/%s.log",  Sim.outexpno);
    fp = fopen(Sim.logfile, "w");
    if (fp == NULL)
      Error("Results directory does not exist");
    else
      fclose(fp);
  }
  
  /* ----------------------------------- initialise user defined global parameters */
  GlobalParameters();
  
  /* -------------------------------------------------------- initialise variables */
  Sim.V.variables = 0;
  Sim.V.scaletime = 1;
  Sim.V.W = NULL;
  Variables(&Sim.V);
  if (Sim.V.variables) {
    Sim.V.W = (struct variable*) malloc(Sim.V.variables*sizeof(struct variable));
    Variables(&Sim.V);
  }

  /* ------------------------------------------------------ initialise time domain */
  ReadParam("int", &Sim.notime, "no_time");
  if (!Sim.notime) {
    ReadParam("char*", timeend,   "time_end");
    ReadParam("char*", timestart, "time_start");
    if (strcmp(timestart, "NA") == 0 || strcmp(timestart, "na") == 0)
      Sim.timestart = -INFINITY;
    else
      Sim.timestart = atof(timestart);
    if (strcmp(timeend, "NA") == 0 || strcmp(timeend, "na") == 0)
      Sim.timeend = INFINITY;
    else
      Sim.timeend = atof(timeend);
  }

  /* -----------------------------find which treatments and individuals to analyse */
  Sim.trtused = integervector(MAXTRTSREPS);
  Sim.indused = integermatrix(MAXTRTSREPS, MAXTRTSREPS);

  ReadParam("char*", ignoreinds,  "ignore");
  ReadParam("char*", treatments,  "treatments");
  ReadParam("char*", individuals, "individuals");

  if (strlen(treatments)) 
    Used(Sim.trtused, treatments);
  else
    for (nt = 0; nt < MAXTRTSREPS; nt++)
      Sim.trtused[nt] = TRUE;

  if (strlen(individuals)) {
    for (nt = 0; nt < MAXTRTSREPS; nt++) if (Sim.trtused[nt])
      Used(Sim.indused[nt], individuals);
  }
  else
    for (nt = 0; nt < MAXTRTSREPS; nt++) if (Sim.trtused[nt])
  for (nr = 0; nr < MAXTRTSREPS; nr++)
    Sim.indused[nt][nr] = TRUE;

  if (strlen(ignoreinds)) Ignore(Sim.indused, ignoreinds);

  /* ------------------------------------------------------------ read in the data */
  ReadParam("char*", Data.filename, "data_file");
  ReadData(&Sim, &Data);

  /* ------------------------------------- initialise user defined McMC parameters */
  Sim.tempfile = charvector(1000);
  Sim.parsfile = charvector(1000);

  ReadParam("int",    &Sim.sim,             "simulate"       );
  ReadParam("int",    &Sim.seed,            "seed"           );
  ReadParam("int",    &Sim.marg,            "marginals"      );
  ReadParam("int",    &Sim.popd,            "prediction"     );
  ReadParam("int",    &Sim.thin,            "thin"           );
  ReadParam("int",    &Sim.range,           "range"          );
  ReadParam("int",    &Sim.samples,         "samples"        );
  ReadParam("int",    &Sim.threads,         "threads"        );
  ReadParam("int",    &Sim.copyind,         "copyind"        );
  ReadParam("int",    &Sim.discard,         "discard"        );
  ReadParam("int",    &Sim.gsl_error,       "gsl_error"      );
  ReadParam("int",    &Sim.population,      "population"     );
  ReadParam("int",    &Sim.covsamples,      "covsamples"     );
  ReadParam("int",    &Sim.phase1adjust,    "phase1adjust"   );
  ReadParam("int",    &Sim.phase2adjust,    "phase2adjust"   );
  ReadParam("int",    &Sim.parsperblock,    "parsperblock"   );
  ReadParam("int",    &Sim.mincovsamples,   "mincovsamples"  );
  ReadParam("int",    &Sim.ignorewarning,   "ignorewarning"  );
  ReadParam("int",    &Sim.usehypermeans,   "usehypermeans"  );
  ReadParam("char*",  Sim.parsfile,         "pars_file"      );
  ReadParam("char*",  Sim.tempfile,         "temp_file"      );
  ReadParam("double", &Sim.xrate,           "xrate"          );
  ReadParam("double", &Sim.lowrate,         "lowrate"        );
  ReadParam("double", &Sim.highrate,        "highrate"       );
  ReadParam("double", &Sim.minprior,        "minprior"       );
  ReadParam("double", &Sim.covthresh,       "covthresh"      );
  ReadParam("double", &Sim.initscale,       "initscale"      );
  ReadParam("double", &Sim.scalefactor,     "scalefactor"    );
  ReadParam("double", &Sim.saturatedlogL,   "saturatedlogL"  );
  ReadParam("double", &Sim.minexchangerate, "minexchangerate");

  /* ----------------------- set gsl error off by default, uncomment to switch on */
  gsl_set_error_handler_off();
  
  if (Sim.threads == 0) Error("threads should equal 1 or more");

  if (Sim.rdp) {
    Sim.mcmc = FALSE;
    Sim.marg = TRUE;
    Sim.popd = FALSE;
    printf("Recalculating derived parameters. McMC is off\n");
  }

  if (!Sim.seed) Sim.seed = time(&t)%100000;

  Sim.discard *= Sim.thin;
  
  /* --------------------------------------------------- set the number of threads */
  omp_set_num_threads(Sim.threads);

  /* printf("%d\n", omp_get_num_devices()); */
  /* exit(0); */
  /* --------------- run on a single processor if the processor argument is passed */
/*   if (argc > 2) { */
/*     int nproc = omp_get_num_procs(); */
/*     int ncores = nproc/2; */

/*     /\* ------------------------------- check that the number of threads is correct *\/ */
/*     if (Sim.threads > ncores) { */
/*       fprintf(stderr, "To run on a single processor reduce threads to %d\n", ncores); */
/*       Error(""); */
/*     } */

/*     int package = atoi(argv[2]); */
/*     if (package != 0 && package != 1) Error("Processor 0 or 1"); */
/*     package *= ncores/2; */

/*     /\* ------------ set affinity mask to force OpenMP to run on a single processor *\/ */
/* #   pragma omp parallel firstprivate(nproc, ncores, package) */
/*     { */
/*       int i; */
/*       int tnum = omp_get_thread_num(); */
/*       kmp_affinity_mask_t mask; */

/*       kmp_create_affinity_mask(&mask); */
/*       for (i = package+tnum/2; i < nproc; i += ncores) */
/*   kmp_set_affinity_mask_proc(i, &mask); */
/*       if (kmp_set_affinity(&mask) != 0) Error("kmp_affinity"); */
/*     } */
/*   } */

  /* -------------------------------------------------------------- begin analysis */
  Bayes(&Sim, &Data);
  RecalculateDerivedParameters(&Sim, &Data);
  Marginals(&Sim, &Data);
  PosteriorPrediction(&Sim, &Data);
  PostAnalysis(&Sim, &Data);
  SimulatedData(&Sim, &Data);
  return SUCCESS;
}

void ReadData(Simulation *Sim, struct data *Data)
{
  int j, k, l, ntrt, nind, np, nt, nr, value, firsttime;
  size_t n;
  double *O, *Oo, x;
  struct data_treatment *T = NULL;
  struct data_individual *R = NULL;
  struct variables *V = &Sim->V;
  char a[1000], *lineptr = NULL;
  FILE *fp;

  fp = fopen(Data->filename, "r");
  if (!fp) Error("no data file");
  Data->ntrt = 0;
  Data->Treatment = NULL;
  Data->maxnind = 0;
  Data->maxnp = 0;
  Data->nvar = -1;

  if (Sim->verbose) printf("Reading in data from %s\n", Data->filename);

  /* ---------------------------------------------------------------- read in data */
  while (fscanf(fp, "%s", a) != EOF)
    if (strncmp(a, "Variables:", 1) == 0) {
      // number of variables
      fscanf(fp, "%d", &Data->nvar);
      Data->nvar++;
    }
    else if (strncmp(a, "Factors:", 1) == 0) {
      // factors in experiment
      if (a[strlen(a)-1] != '\n') {
  n = getline(&lineptr, &n, fp);
  lineptr[n-1] = '\0';
  strcpy(Sim->factors, lineptr);
  free(lineptr);
    lineptr = NULL;
      }
    }
    else if (strncmp(a, "Treatment", 1) == 0) {
      // new treatment group
      if (Data->nvar == -1) Error("Set Variables first in data file");

      // increment no. of treatments if first occurence or last treatment group
      // contained individuals
      if (T == NULL || T->nind > 0) ntrt = ++Data->ntrt;
      Data->Treatment = (struct data_treatment*) 
  realloc(Data->Treatment, ntrt*sizeof(struct data_treatment));
      T = &Data->Treatment[ntrt-1];
      T->nind = 0;
      T->usevar = integervector(Data->nvar);
      for (j = 1; j < Data->nvar; j++)
  fscanf(fp, "%d", &T->usevar[j]);
      T->Individual = NULL;
      T->levels[0] = '\0';
    }
    else if (strncmp(a, "Levels:", 1) == 0) {
      // levels for this treatment
      if (T == NULL) Error("Data file in wrong format");
      if (a[strlen(a)-1] != '\n') {
  n = getline(&lineptr, &n, fp);
  lineptr[n-1] = '\0';
  strcpy(T->levels, lineptr);
  free(lineptr);
    lineptr = NULL;
      }
    }
    else if (strncmp(a, "Replicate", 1) == 0 || strncmp(a, "Individual", 1) == 0) {
      // new individual
      if (Data->nvar == -1) Error("Set Variables first in data file");
      if (T == NULL) Error("Data file in wrong format");
      nind = ++T->nind;
      if (Data->maxnind < nind) Data->maxnind = nind;
      T->Individual = (struct data_individual*)
  realloc(T->Individual, nind*sizeof(struct data_individual));
      R = &T->Individual[nind-1]; 
      R->n = 0;
      R->np = 0;
      R->Y = doublematrix(1, Data->nvar);
      R->value = uintegermatrix(1, Data->nvar);
    }
    else if (a[0] == '%') {
      // a comment, if last character of a is not a new line, read rest of line
      if (a[strlen(a)-1] != '\n') {
  getline(&lineptr, &n, fp);
  free(lineptr);
  lineptr = NULL;
      }
    }
    else {
      if (strncmp(a, "#", 1) == 0) {
  // a commented out time point for all variables
  value = NA;
  fscanf(fp, "%s", a);
      }
      else
  // a new time point
  value = MEASURED;

      if (R == NULL) Error("Data file in wrong format");
      np = R->np;
      R->Y[np][0] = x = atof(a);
      for (j = 1; j < Data->nvar; j++) {
  fscanf(fp, "%s", a);
  if (T->usevar[j]) {
    if (strcmp(a, "NA") == 0 || strcmp(a, "na") == 0 ||
        strcmp(a, "ND") == 0 || strcmp(a, "nd") == 0) {
      R->value[np][j] = NA;
      R->Y[np][j] = NAN;
    }
    else if (strcmp(a, "BD") == 0 || strcmp(a, "bd") == 0) {
      R->value[np][j] = BELOW_DETECTION;
      R->Y[np][j] = 0;
      if (T->usevar[j] == 1) R->n++;
    }
    else if (strcmp(a, "AD") == 0 || strcmp(a, "ad") == 0) {
      R->value[np][j] = ABOVE_DETECTION;
      R->Y[np][j] = 0;
      if (T->usevar[j] == 1) R->n++;
    }
    else {
      R->value[np][j] = value;
      R->Y[np][j] = atof(a);
      if (T->usevar[j] == 1) R->n++;
    }
  }
  else
    R->value[np][j] = NA;
      }

      if (Sim->notime || (Sim->timestart <= x && x <= Sim->timeend)) {
  np = ++R->np;
  if (np > Data->maxnp) Data->maxnp = np;
  R->Y = (double**) realloc(R->Y, (np+1)*sizeof(double*));
  R->Y[np] = (double*) calloc(Data->nvar, sizeof(double));
  R->value = (uint**) realloc(R->value, (np+1)*sizeof(uint*));
  R->value[np] = (uint*) calloc(Data->nvar, sizeof(uint));
      }
    }
  fclose(fp);

  // remove last treatment if it contains no individuals
  if (T->nind == 0) Data->ntrt--;

  /* ----------------------------------------- find range of treatments to analyse */
  nt2 = MAXTRTSREPS-1;
  while(!Sim->trtused[nt2] || nt2 >= Data->ntrt) nt2--;
  nt1 = 0;
  while(!Sim->trtused[nt1] && nt1 <= nt2) nt1++;
  if (nt1 > nt2) Error("miss-specificiation of treatments to analyse");

  /* ------------------------------------------------------- output converted data */
  sprintf(a, "Results/%s.data", Sim->outexpno);
  fp = fopen(a, "w");
  if (V->variables) {
    O = doublevector(V->variables);
    Oo = doublevector(V->variables);
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++)if (Sim->indused[nt][nr]) {
        fprintf(fp, "data %d %d\n", nt, nr);
        R = &Data->Treatment[nt].Individual[nr];
        for (l = 1; l < V->variables; l++) {
          firsttime = TRUE;
          for (k = 0; k < R->np; k++) {
            SetNAN(O, V->variables);
            SetNAN(Oo, V->variables);
            OutputData(nt, nr, O, Oo, R->Y[k], R->value[k], k);
            if (!isnan(O[l])) {
              if (firsttime) {
                fprintf(fp, "%d\n", l);
                firsttime = FALSE;
              }
              fprintf(fp, "%f %.3e %.3e\n", R->Y[k][0]*V->scaletime, O[l], Oo[l]);
            }
          }
        }
      }
    free(O); free(Oo);
  }
  else {
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
      T = &Data->Treatment[nt];
      for (nr = nr1; nr <= NR2; nr++)if (Sim->indused[nt][nr]) {
        fprintf(fp, "data %d %d\n", nt, nr);
  R = &T->Individual[nr];
  for (k = 0; k < R->np; k++) {
    for (l = 0; l < Data->nvar; l++)
      if (T->usevar[l] || l == 0) fprintf(fp, "%.1e ", R->Y[k][l]);
    fprintf(fp, "\n");
  }
      }
    }
  }
  fclose(fp);
}

void Var(struct variables *V, int index, char *name) 
{
  struct variable *W = V->W;
  
  if (W == NULL) {
    V->variables = max(V->variables, index+1);
    return;
  }
  strcpy(W[index].name, name);
}  

void ScaleTime(struct variables *V, double scale)
{
  V->scaletime = scale;
}

void Create(Hyperparameters *H, int *h, int block, int par,
      Matrix M, char *name, int type, double p1, double p2)
{
  int i, j;
  CreateHyperBlock(H, block);
  for (i = 0; i < M.n; i++) {
    Hyper(H, *h, Real, name, type, p1, p2);
    for (j = 0; j < M.m; j++)
      AddHyperToATreatment(H, block, *h, par, M.a[i][j]);
    (*h)++;
  }
}

void Hyper(Hyperparameters *HP, int index, int number, char *name,
     int distribution, double shape_or_loc, double scale)
{
  if (index < 0) Error("hyper parameter index must be >= 0");

  /* ---------------------- first pass init=TRUE, count number of hyper parameters */
  if (HP->init) {
    HP->nhyper_parameters = max(HP->nhyper_parameters, index+1);
    return;
  }

  /* ----------------------------------- second pass init=FALSE, initialise hypers */
  if (number != REAL && number != INTEGER) {
    fprintf(stderr, "hyper parameter %d does not have a valid number\n", index);
    Error("");
  }
  if (distribution != NOTAPPLICABLE &&
      distribution != NORMAL        && distribution != GAMMA    &&
      distribution != EXPONENTIAL   && distribution != BINOMIAL &&
      distribution != POISSON       && distribution != SICHI2   &&
      distribution != INVGAMMA      && distribution != PARETO   &&
      distribution != BETA) {
    fprintf(stderr, "hyper parameter %d does not have a valid distribution\n", 
      index);
    Error("");
  }
  if (distribution != NOTAPPLICABLE && scale <= 0) {
    fprintf(stderr, "hyper parameter %d must have scale > 0\n", index);
    Error("");
  }
  if (HP->H[index].index != -1) {
    fprintf(stderr, "hyper parameter %d already used\n", index);
    Error("");
  }

  strcpy(HP->H[index].name, name);
  HP->H[index].index = index;
  HP->H[index].number = number;
  HP->H[index].hdist = distribution;
  HP->H[index].shape = shape_or_loc;
  HP->H[index].location = shape_or_loc;
  HP->H[index].scale = scale;
}

void CreateHyperBlock(Hyperparameters *HP, int index)
{
  // nhypers_used = number of hypers associated with this block
  // hyper_used = list of hypers associated with this block
  // nassociated_pars = 2
  // associated_par = {14, 16}
  // ntreatments_with_parameter = {16, 8}
  // treatment_with_parameter[0] = {0...15}
  // treatment_with_parameter[1] = {8...15}
  // nhypers_with_par_treatment[0] = {4, 7, 3, ... # of treatments}
  // nhypers_with_par_treatment[1] = {4, 7, 3, ... # of treatments}
  // hyper_with_par_treatment[0][0] = {0, 3, 4, 5}
  // hyper_with_par_treatment[0][#trts] = {1, 6, 8}
  // hyper_with_par_treatment[1][0] = {0, 3, 4, 5}
  // hyper_with_par_treatment[1][#trts] = {1, 6, 8}
  // dependent_block[0] = 2
  // dependent_block[1] = 3

  if (index < 0) Error("hyper set index must be >= 0");

  if (HP->init) {
    HP->nblocks = max(HP->nblocks, index+1);
    return;
  }

  if (HP->block[index].use == TRUE) {
    fprintf(stderr, "Hyper block %d created twice\n", index);
    Error("");
  }

  HP->block[index].use = TRUE;  
  HP->block[index].nhypers_used = 0;
  HP->block[index].nassociated_pars = 0;
  HP->block[index].associated_par = NULL;
  HP->block[index].ntreatments_with_parameter = NULL;
  HP->block[index].treatment_with_parameter = NULL;
  HP->block[index].hdistribution = NOTAPPLICABLE;
  HP->block[index].pdistribution = NOTAPPLICABLE;

  HP->block[index].hyper_used = integervector(HP->nhyper_parameters);
  HP->block[index].nhypers_with_treatment = NULL;
  HP->block[index].hyper_with_treatment = NULL;
}

void AddHyperToATreatment(Hyperparameters *HP, int index, int h, int pindex, int nt)
{
  if (HP->init) return;

  int i, j;
  Hyperblock *HB = &HP->block[index];

  if (nt > nt2 || !HP->trtused[nt]) {
    fprintf(stderr, "AddHyperToATreatment: hyper %d not added to treatment %d as "
      "this treatment is not used\n", h, nt);
    return;
  }

  /* ------------- check that distribution of hyper h matches others in same block */
  if (HB->hdistribution == NOTAPPLICABLE) 
    HB->hdistribution = HP->H[h].hdist;
  /* else if (HB->hdistribution != HP->H[h].hdist) { */
  /*   fprintf(stderr, "distribution of hyper parameter %d does not match " */
  /*       "with others in hyper block %d", h, index); */
  /*   Error(""); */
  /* } */

  /* ----------------------------------- add hyper to block if not already done so */
  if (In_set(h, HB->hyper_used, HB->nhypers_used) == -1)
    HB->hyper_used[HB->nhypers_used++] = h;

  /* --------------------------------------- associate parameter pindex with block */
  i = In_set(pindex, HB->associated_par, HB->nassociated_pars);
  /*- ---------------------------------------- parameter pindex not yet associated */
  if (i == -1) {
    i = ++HB->nassociated_pars;
    HB->associated_par = (int*) realloc(HB->associated_par, i*sizeof(int));
    HB->ntreatments_with_parameter = (int*) realloc(HB->ntreatments_with_parameter, 
               i*sizeof(int));
    HB->treatment_with_parameter = (int**) realloc(HB->treatment_with_parameter, 
               i*sizeof(int*));
    HB->nhypers_with_treatment = (int**) realloc(HB->nhypers_with_treatment,
             i*sizeof(int*));
    HB->hyper_with_treatment = (int***) realloc(HB->hyper_with_treatment,
            i*sizeof(int**));
    i--;
    HB->associated_par[i] = pindex;
    HB->ntreatments_with_parameter[i] = 0;
    HB->treatment_with_parameter[i] = NULL;
    HB->nhypers_with_treatment[i] = integervector(nt2+1);
    HB->hyper_with_treatment[i] = integermatrix(nt2+1, HP->nhyper_parameters);
  }
  
  /* ------------- add hyper h to (parameter i, treatment nt) combination in block */
  if (In_set(h, HB->hyper_with_treatment[i][nt], HB->nhypers_with_treatment[i][nt])
      != -1) {
    fprintf(stderr, "treatment %d already added to hyperparameter %d for "
      "hyper block %d\n", nt, h, index);
    Error("");
  }
  HB->hyper_with_treatment[i][nt][HB->nhypers_with_treatment[i][nt]++] = h;

  /* ------------------------- add treatment nt to associated parameter k in block */
  j = In_set(nt, HB->treatment_with_parameter[i], HB->ntreatments_with_parameter[i]);
  if (j == -1) {
    j = ++HB->ntreatments_with_parameter[i];
    HB->treatment_with_parameter[i] = (int*) realloc(HB->treatment_with_parameter[i], 
                 j*sizeof(int));
    j--;
  }
  HB->treatment_with_parameter[i][j] = nt;
}

void AddHyperToTreatments(Hyperparameters *HP, int index, int h, int pindex, ...)
{
  if (HP->init) return;

  int i, j, nt;
  Hyperblock *HB = &HP->block[index];
  va_list ap;

  /* --------------------------------------- associate parameter pindex with block */
  i = In_set(pindex, HB->associated_par, HB->nassociated_pars);
  /* ----------------------------------------- parameter pindex not yet associated */
  if (i == -1) {
    i = ++HB->nassociated_pars;
    HB->associated_par = (int*) realloc(HB->associated_par, i*sizeof(int));
    HB->ntreatments_with_parameter = (int*) realloc(HB->ntreatments_with_parameter, 
               i*sizeof(int));
    HB->treatment_with_parameter = (int**) realloc(HB->treatment_with_parameter, 
               i*sizeof(int*));
    HB->nhypers_with_treatment = (int**) realloc(HB->nhypers_with_treatment,
             i*sizeof(int*));
    HB->hyper_with_treatment = (int***) realloc(HB->hyper_with_treatment,
            i*sizeof(int**));
    i--;
    HB->associated_par[i] = pindex;
    HB->ntreatments_with_parameter[i] = 0;
    HB->treatment_with_parameter[i] = NULL;
    HB->nhypers_with_treatment[i] = integervector(nt2+1);
    HB->hyper_with_treatment[i] = integermatrix(nt2+1, HP->nhyper_parameters);
  }
  
  /* ------------- check that distribution of hyper h matches others in same block */
  if (HB->hdistribution == NOTAPPLICABLE)
    HB->hdistribution = HP->H[h].hdist;
  /* else if (HB->hdistribution != HP->H[h].hdist) { */
  /*   fprintf(stderr, "distribution of hyper parameter %d does not match " */
  /*       "with others in hyper block %d", h, index); */
  /*   Error(""); */
  /* } */

  /* ------------------------------------------- add hypers to treatments in block */
  va_start(ap, pindex);
  while ((nt = va_arg(ap, int)) != EOL) {
    if (nt > nt2 || !HP->trtused[nt]) {
      fprintf(stderr, "AddHyperToATreatments: hyper %d not added to treatment %d "
          "as this treatment is not used\n", h, nt);
      continue;
    }

    /* --------------------------------- add hyper to block if not already done so */
    if (In_set(h, HB->hyper_used, HB->nhypers_used) == -1)
      HB->hyper_used[HB->nhypers_used++] = h;
  

    /* ----------- add hyper h to (parameter i, treatment nt) combination in block */
    if (In_set(h, HB->hyper_with_treatment[i][nt], HB->nhypers_with_treatment[i][nt])
  != -1) {
      fprintf(stderr, "treatment %d already added to hyperparameter %d for "
        "hyper block %d\n", nt, h, index);
      Error("");
    }
    HB->hyper_with_treatment[i][nt][HB->nhypers_with_treatment[i][nt]++] = h;

    /* ----------------------- add treatment nt to associated parameter k in block */
    j = In_set(nt, HB->treatment_with_parameter[i],
         HB->ntreatments_with_parameter[i]);
    if (j == -1) {
      j = ++HB->ntreatments_with_parameter[i];
      HB->treatment_with_parameter[i] =
  (int*) realloc(HB->treatment_with_parameter[i], j*sizeof(int));
      j--;
    }
    HB->treatment_with_parameter[i][j] = nt;
  }
  va_end(ap);

}

void Block(struct parameterset *Q, int npar, ...)
{
  int block, i, k, start, end, *pars;
  struct blocks *B = &Q->Blocks;
  va_list ap;

  va_start(ap, npar);
  if (npar == 0) {
    start = va_arg(ap, int);
    end = va_arg(ap, int);
    npar = end-start+1;
    pars = integervector(npar);
    k = 0;
    for (i = start; i <= end; i++)
      pars[k++] = i;
  }
  else {
    pars = integervector(npar);
    for (i = 0; i < npar; i++)
      pars[i] = va_arg(ap, int);
  }
  va_end(ap);

  /* ------------------------------------ if a parameter is already set do nothing */
  for (block = 0; block < B->nblocks; block++)
    for (i = 0; i < npar; i++)
      if (pars[i] == B->pblock[block][1]) {
  free(pars);
  return;
      }
  
  B->pblock = (int**) realloc(B->pblock, (B->nblocks+1)*sizeof(int*));
  B->pblock[B->nblocks] = integervector(npar+1);
  B->pblock[B->nblocks][0] = npar;

  for (i = 1; i <= npar; i++)
    B->pblock[B->nblocks][i] = pars[i-1];

  B->nblocks++;
  free(pars);
}

void Par(struct parameterset *Q, int index, int number, char* name, int type, ...) 
{
  /* ---------------------------------------- need go no further after second call */
  if (Q->pass == 3) return;

  /* ---------------------- if second pass and parameter already set through error */
  struct parameter *P = Q->P;
  if (Q->pass == 2 && P[index].type) {
    fprintf(stderr, "Parameter %d already set", index);
    Error("");
  }

  /* ------------------------------------------------ first call, count parameters */
  if (P == NULL) {
    Q->parameters = max(Q->parameters, index+1);
    if (type != CONSTANT && type != INDIVIDUAL_CONSTANT) Q->npar++;
    return;
  }

  int i, distribution;
  double sum;
  va_list ap;

  if (index >= MAXPARS) {
    fprintf(stderr, "Increase MAXPARS in bayes.h");
    Error("");
  }

  /* --------------------------------- second call, initialise prior distributions */
  if (number != REAL && number != INTEGER && number != CATEGORY) {
    fprintf(stderr, "parameter %d does not have a valid number\n", index);
    Error("");
  }
  if (type != NORMAL 
      && type != GAMMA 
      && type != SICHI2 
      && type != UNIFORM 
      && type != CIRCULAR
      && type != PARETO 
      && type != TRUNCNORMAL 
      && type != EXPONENTIAL 
      && type != BINOMIAL 
      && type != POISSON 
      && type != INVGAMMA 
      && type != UNIFORM0 
      && type != BETA 
      && type != NORMALMEAN
      && type != CATEGORICAL 
      && type != GEOMETRIC 
      && type != USERDEFINED 
      && type != HYPER 
      && type != CONSTANT 
      && type != INDIVIDUAL_CONSTANT 
      && type != AUXILLARY 
      && type != AUXILLARY_HYPER 
      && type != DERIVED) {
    fprintf(stderr, "parameter %d does not have a valid type or prior\n", index);
    Error("");
  }

  P[index].index = index;
  P[index].number = number;
  strcpy(P[index].name, name);

  if (type == AUXILLARY_HYPER) {
    P[index].type = AUXILLARY;
    P[index].hashyper = TRUE;
    P[index].distribution = NORMAL;
    P[index].circular = FALSE;
  }
  else if (type == CONSTANT || type == INDIVIDUAL_CONSTANT || type == DERIVED || 
     type == AUXILLARY) {
    P[index].type = type;
    P[index].hashyper = FALSE;
    P[index].distribution = NOTUSED;
    P[index].circular = FALSE;
  }
  else {
    P[index].type = FITTED;
    if (type == HYPER) {
      P[index].hashyper = TRUE;
      if (number == INTEGER) {
  fprintf(stderr, "parameter %d has a hyper so it must be a Real\n", index);
  Error("");
      }
    }
    else
      P[index].hashyper = FALSE;
    
    if (type == CIRCULAR) {
      P[index].circular = TRUE;
      type = UNIFORM;
    }
    else
      P[index].circular = FALSE;

  }
  
  // nothing more to do if the parameter is derived, auxillary
  if (type == DERIVED || type == AUXILLARY || type == INDIVIDUAL_CONSTANT) return;

  va_start(ap, type);
  if (type == CONSTANT) {
    if (number == REAL)
      P[index].location = va_arg(ap, double);
    else
      P[index].location = va_arg(ap, int);
    /* if (number == REAL && P[index].location == 0) { */
    /*   fprintf(stderr, "WARNING use decimal points in Parameters() for Reals: " */
    /*         "parameter %d = 0\n", index); */
    /*   Error(""); */
    /* } */
    va_end(ap);
    return;
  }

  // should only be here if the parameter has a defined prior distribution
  if (P[index].hashyper)
    P[index].distribution = distribution = va_arg(ap, int); 
  else
    P[index].distribution = distribution = type;

  if (distribution == EXPONENTIAL) {
    if (number != REAL) {
      fprintf(stderr, "parameter %d: must be a Real", index);
      Error("");
    }
    if (type == HYPER)
      P[index].hyper_location_block_index = va_arg(ap, int);
    else {
      P[index].location = va_arg(ap, double);
      P[index].scale = P[index].location*P[index].location;
    }
  }
  else if (distribution == POISSON) {
    if (number != INTEGER) {
      fprintf(stderr, "parameter %d: must be a Integer", index);
      Error("");
    }
    if (type == HYPER)
      P[index].hyper_location_block_index = va_arg(ap, int);
    else {
      P[index].location = va_arg(ap, double);
      P[index].scale = P[index].location;
    }
  }
  else if (distribution == GEOMETRIC) {
    if (number != INTEGER) {
      fprintf(stderr, "parameter %d: must be a Integer", index);
      Error("");
    }
    if (type == HYPER)
      P[index].hyper_location_block_index = va_arg(ap, int);
    else {
      P[index].location = va_arg(ap, double);
      P[index].scale = P[index].location;
      if (P[index].location < 0 || P[index].location > 1) {
  fprintf(stderr, "parameter %d: value must be between 0 and 1", index);
  Error("");
      }
    }
  }
  else if (distribution == UNIFORM0)
    if (type == HYPER)
      P[index].hyper_scale_block_index = va_arg(ap, int);
    else
      if (number == INTEGER)
  P[index].scale = va_arg(ap, int);
      else
  P[index].scale = va_arg(ap, double);
  else if (distribution == GAMMA || distribution == INVGAMMA) {
    if (number != REAL) {
      fprintf(stderr, "parameter %d: must be a Real", index);
      Error("");
    }
    if (type == HYPER) {
      P[index].shape = va_arg(ap, double);
      P[index].hyper_scale_block_index = va_arg(ap, int);
    }
    else {
      P[index].shape = va_arg(ap, double);
      P[index].scale = va_arg(ap, double);
    }
  }
  else if (distribution == BINOMIAL) {
    if (number != INTEGER) {
      fprintf(stderr, "parameter %d: must be an Integer", index);
      Error("");
    }
    if (type == HYPER) {
      P[index].shape = va_arg(ap, int);
      P[index].hyper_scale_block_index = va_arg(ap, int);
    }
    else {
      P[index].shape = va_arg(ap, int);
      P[index].scale = va_arg(ap, double);
      if (P[index].scale < 0 || P[index].scale > 1) {
  fprintf(stderr, "parameter %d must have 0 <= p <= 1\n", index);
  Error("");
      }
    }
    if (P[index].shape < 0) {
      fprintf(stderr, "parameter %d must have n > 0\n", index);
      Error("");
    }
  }
  else if (distribution == NORMALMEAN) {
    if (type != HYPER) {
      fprintf(stderr, "parameter %d: must have the mean as a "
        "hyperparameter\n", index);
      Error("");
    }
    P[index].hyper_location_block_index = va_arg(ap, int);
    P[index].scale = va_arg(ap, double);
    if (P[index].scale <= 0) {
      fprintf(stderr, "parameter %d must have scale > 0\n", index);
      Error("");
    }
  }
  else if (distribution == CATEGORICAL) {
    if (number != CATEGORY) {
      fprintf(stderr, "parameter %d: must be a Category", index);
      Error("");
    }
    P[index].ncategories = va_arg(ap, int);
    P[index].p = doublevector(P[index].ncategories);
    sum = 0;
    for (i = 0; i < P[index].ncategories-1; i++) {
      P[index].p[i] = va_arg(ap, double);
      if (P[index].p[i] < 0) {
  fprintf(stderr, "parameter %d: Category %d is negative", index, i);
  Error("");
      }
      sum += P[index].p[i];
    }
    if (sum > 1.) {
      fprintf(stderr, "parameter %d: Categories sum to > 1", index);
      Error("");
    }
  }
  else if (distribution == BETA) {
    if (number != REAL) {
      fprintf(stderr, "parameter %d: must be a Real", index);
      Error("");
    }
    P[index].shape = va_arg(ap, double);
    P[index].scale = va_arg(ap, double);
  }
  else if (distribution != USERDEFINED) {// eg NORMAL
    if (type == HYPER || type == AUXILLARY_HYPER) {
      P[index].hyper_location_block_index = va_arg(ap, int);
      P[index].hyper_scale_block_index = va_arg(ap, int);
    }
    else if (number == REAL) {
      P[index].location = va_arg(ap, double);
      P[index].scale = va_arg(ap, double);
    }
    else {
      P[index].location = va_arg(ap, int);
      P[index].scale = va_arg(ap, int);
    }
  }
  va_end(ap);
}  

void Bayes(Simulation *Sim, struct data *Data)
{
  int r, c, i, j, k, nt, nr, flag, **redotemp, nadj, memory;
  int np, block, *nonblocked, *pb, npars;
  int block1, block2;
  int v = Sim->verbose, newblock;
  double f, df, *lnL, **y;
  double logML, seML, mintemp, minX, X, tmpt[MAXTRTSREPS];
  double tmptheta[MAXPARS], *mean, *var, *hmean, *error2;
  char filename[100], backup[1000];
  FILE *fp, *fpl;
  struct sim_treatment *T;
  struct sim_individual *Ind;
  struct data_individual *Dind;
  struct chain *Chain, *Chain0;
  struct variables *V = &Sim->V;
  struct parameterset *Q;
  struct parameter_block *PB;
  struct blocks *B;
  struct line ElogL;
  struct hyper *H;
  Hyperparameters *HP = &Sim->HP;
  Hyperblock *HB;
  TimePoints TP;
  double temperatures[NTEMP] = {
    1, 0.65, 0.45, 0.30, 0.20, 0.13, 0.09, 0.062, 0.046, 0.036, 0.024, 0.015, 
    0.008, 0.0048, 0.0030, 0.0019, 0.0010, 0.0006, 0.00035, 0.00020, 0.00010, 
    0.00004, 0.000016, 0.000006, 0.0000006, 6e-30, 6e-100, 6e-200, 6e-300};
  
  Sim->maxnp = 0;
  
  if (Data->ntrt == 0) {
    fprintf(stderr, "Data file contains no treatments\n");
    Error("");
  }
  Sim->Treatment = (struct sim_treatment*) 
    malloc(Data->ntrt*sizeof(struct sim_treatment));
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    if (Data->Treatment[nt].nind == 0) {
      fprintf(stderr, "treatment %d contains no members. "
        "Remove this treatment from the data file\n", nt);
      Error("");
    }
    Sim->Treatment[nt].Individual = (struct sim_individual*) 
      malloc(Data->Treatment[nt].nind*sizeof(struct sim_individual));
    }

  memory = Data->ntrt*sizeof(struct sim_treatment);
  if (v == 2) printf("m1: %d\n", memory);

  if (v == 2) printf("init hypers\n");
  /* --------------- initialise hyperparameters and associate them with treatments */
  /* ---- first pass: find values of Sim->HP.nhyper_parameters and Sim->HP.nblocks */
  HP->nhyper_parameters = 0;
  HP->nblocks = 0;
  HP->init = TRUE;

  Treatments Tr;
  Tr.n = 0;
  Tr.treatment = integervector(nt2+1);
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      Tr.treatment[Tr.n++] = nt;
  HyperParameters(Tr, &Sim->HP);

  if (HP->nhyper_parameters) {
    HP->H = (struct hyper*) malloc(HP->nhyper_parameters*sizeof(struct hyper));
    for (i = 0; i < HP->nhyper_parameters; i++) {
      HP->H[i].use = FALSE;
      HP->H[i].index = -1;
    }

    if (HP->nblocks) {
      HP->block = (Hyperblock*) malloc(HP->nblocks*sizeof(Hyperblock));
      for (i = 0; i < HP->nblocks; i++) {
  HP->block[i].use = FALSE;
  HP->block[i].nhypers_used = 0;
      }
    }

    HP->nind = integervector(nt2+1);
    /* ------------------------- count the number of individuals in each treatment */
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
  HP->nind[nt] = 0;
  for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
    HP->nind[nt]++;
    }

    /* -------------------------------------------- second pass: set up everything */
    HP->init = FALSE;
    HP->ntrt = Data->ntrt;
    HP->trtused = Sim->trtused;
    HyperParameters(Tr, &Sim->HP);
  }
  H = HP->H;

  if (v == 2) printf("init ind pars\n");
  /* -------------------------------------------- initialise individual parameters */
  int maxnp = 0;
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    T = &Sim->Treatment[nt];
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      Q = &Ind->Q;
      Q->P = NULL;
      Q->npar = 0;
      Q->parameters = 0;
      Q->nauxillary = 0;
      Q->ninteger = 0;
      B = &Q->Blocks;
      B->nblocks = 0;
      B->pblock = NULL;

      /* ----------------------------- count parameters and declare initial memory */
      Q->pass = 1;
      Parameters(nt, nr, tmptheta, Q);
      np = Q->parameters;

      Q->P = (struct parameter*) malloc(np*sizeof(struct parameter));
      for (k = 0; k < np; k++) {
  Q->P[k].type = NOTUSED;
  Q->P[k].hashyper = FALSE;
      }

      /* ------------------------------------------ initialise prior distributions */
      Q->pass = 2;
      Parameters(nt, nr, tmptheta, Q);
      Q->pass = 3;

      /* ------------------------------------- make list of all integer parameters */
      for (k = 0; k < np; k++) 
  if (Q->P[k].number == INTEGER) Q->ninteger++;
      if (Q->ninteger) {
  Q->integer = integervector(Q->ninteger);
  for (i = 0, k = 0; k < np; k++) 
    if (Q->P[k].number == INTEGER) Q->integer[i++] = k;
      }

      /* ----------------------------------- make list of all auxillary parameters */
      for (k = 0; k < np; k++) 
  if (Use(AUXILLARY, Q->P[k].type)) Q->nauxillary++;
      if (Q->nauxillary) {
  Q->auxillary = integervector(Q->nauxillary);
  for (i = 0, k = 0; k < np; k++) 
    if (Use(AUXILLARY, Q->P[k].type)) Q->auxillary[i++] = k;
      }

      /* ------------------------------------ find individual with most parameters */
      if (np > Sim->maxnp) Sim->maxnp = np;
      
      int cnp = 0;
      for (k = 0; k < np; k++) 
  if (Use(ALL, Q->P[k].type)) cnp++;
      if (cnp > maxnp) {
  maxnp = cnp;
  Sim->ntp = nt;
  Sim->nrp = nr;
      }
    }
  }

  if (v == 2) printf("init hypers\n");

  /* ------------------------------------------------- initialise hyper parameters */
  /* ----------------------------- using the paramter declarations in Parameters() */
  if (HP->nblocks) {
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
      nr = 0;
      while (!Sim->indused[nt][nr]) nr++;
  
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      for (k = 0; k < Q->parameters; k++) 
  if (Q->P[k].hashyper) {
    if (Q->P[k].distribution == GAMMA    || Q->P[k].distribution == INVGAMMA ||
        Q->P[k].distribution == UNIFORM0 || Q->P[k].distribution == BINOMIAL)
      block1 = Q->P[k].hyper_scale_block_index;
    else
      block1 = Q->P[k].hyper_location_block_index;
    SetupHyperBlock(HP, Q, k, nt, block1);
    
    if (Q->P[k].distribution != EXPONENTIAL && 
        Q->P[k].distribution != POISSON     &&
        Q->P[k].distribution != UNIFORM0    &&
        Q->P[k].distribution != NORMALMEAN  &&
        Q->P[k].distribution != GAMMA) {
      block2 = Q->P[k].hyper_scale_block_index;
      SetupHyperBlock(HP, Q, k, nt, block2);
      
      j = In_set(k, HP->block[block1].associated_par, 
           HP->block[block1].nassociated_pars);
      if (j == -1) {
        fprintf(stderr, "parameter %d not in block %d", k, block1);
        Error("");
      }

      HP->block[block1].dependent_block[j] = block2;
      
      j = In_set(k, HP->block[block2].associated_par, 
           HP->block[block2].nassociated_pars);
      if (j == -1) {
        fprintf(stderr, "parameter %d not in block %d", k, block2);
        Error("");
      }
      HP->block[block2].dependent_block[j] = block1;
    }
    //  else
    //    HP->block[block1].dependent_block[j] = -1;
    
  }
    }
    /* ------------------------------------------------ remove unused hyper blocks */
    for (i = 0; i < HP->nblocks; i++) {
      HB = &HP->block[i];
      if (HB->use == 1) HB->use = FALSE;
      
      if (HB->use) {
  /* -------------- sort ascending hyper_used vector: must be done for Gibbs */
  Int_HeapSort(HB->hyper_used-1, HB->nhypers_used);

  /* - and only use hyperparameters that are associated with used treatments */
  for (j = 0; j < HB->nassociated_pars; j++)
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (k = 0; k < HB->nhypers_with_treatment[j][nt]; k++)
        H[HB->hyper_with_treatment[j][nt][k]].use = TRUE;
      }
    }
  }

  if (v == 2) printf("init ladders\n");
  /* ------------------------------------- initialise and read temperature ladders */
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
      Sim->Treatment[nt].Individual[nr].temp = NULL;

  if (Sim->population) {
    flag = FALSE;
    if ((fpl = fopen(Sim->tempfile, "r"))) {
      /* --------------------------------------------- read temperatures from file */
      while ((j = fscanf(fpl, "%d%d%d", &nt, &nr, &i)) != EOF) {
  if (j != 3) {
    fprintf(stderr, "%s is in wrong format: treatment individual "
      "no_of_chains", Sim->tempfile);
    Error("");
  }

  for (c = 0; c < i; c++)
    if (!fscanf(fpl, "%lf", &tmpt[c])) {
      fprintf(stderr, "%s is in wrong format for treatment %d, individual "
        "%d\n", Sim->tempfile, nt, nr);
      Error("");
    }

  if (nt > NT2 || nr > NR2) {
    if (!flag)
      fprintf(stderr, "WARNING: %s does not match %s\n", Sim->tempfile,
        Data->filename);
    flag = TRUE;
  }
  else if (Sim->trtused[nt] && Sim->indused[nt][nr]) {
    /* ------------------------------ add further lower default temperatures */
    i = 0;
    while (i < NTEMP && temperatures[i++] > tmpt[c-1]);
    while (i < NTEMP) tmpt[c++] = temperatures[i++];
    
    /* -------------------------------------------------------- setup ladder */
    Ind = &Sim->Treatment[nt].Individual[nr];
    Ind->chains = c+1;
    Ind->temp = doublevector(Ind->chains);
    for (c = 0; c < Ind->chains-1; c++)
      Ind->temp[c] = tmpt[c];
    Ind->temp[Ind->chains-1] = 0;
  }
      }
      fclose(fpl);
    }

    /* --------------- use default temperatures for individuals not in .temp file */
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  if (Ind->temp == NULL) {
    if (Sim->verbose) fprintf(stderr, "using default temperature ladder "
            "for treatment %d, individual %d\n", nt, nr);
    Ind->chains = NTEMP;
    Ind->temp = doublevector(Ind->chains);
    for (c = 0; c < Ind->chains-1; c++)
      Ind->temp[c] = temperatures[c];
    Ind->temp[Ind->chains-1] = 0;
  }
      }


    /* --------------------Sim->chains is set to the length of the longest ladder */
    Sim->chains = 0;
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  if (Ind->chains > Sim->chains) {
    Sim->chains = Ind->chains;
    i = nt;
    j = nr;
  }
      }
    /* ------------------------ use only one ladder if there are hyper parameters */
    if (HP->nhyper_parameters) {
      /* ------------------------------- copy longest ladder to all other ladders */
      for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
  for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
    if (nt == i && nr == j) continue;
    Ind = &Sim->Treatment[nt].Individual[nr];
    Ind->chains = Sim->chains;
    Ind->temp = (double*) realloc(Ind->temp, Sim->chains*sizeof(double));
    for (c = 0; c < Sim->chains; c++)
      Ind->temp[c] = Sim->Treatment[i].Individual[j].temp[c];
  }
    }  
  }
  else {
    /* ---------------------------- population=0 so no temperature ladder required */
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  Ind->temp = doublevector(1);
  Ind->temp[0] = 1;
  Ind->chains = 1;
      }
  }
    
  if (v == 2) printf("init individual parameters\n");
  /* ------------------------------------------- initialise individual parameters */
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Ind->Chain = (struct chain*) malloc(Ind->chains*sizeof(struct chain));
      Ind->chainno = integervector(MAXSTORE);
      Q = &Ind->Q;
      np = Q->parameters;
      B = &Q->Blocks;

      if (v == 2) printf("init chains\n");
      for (c = 0; c < Ind->chains; c++) {
        Chain = &Ind->Chain[c];
        Chain->logLstore = doublevector(MAXSTORE);
        Chain->logpoststore = doublevector(MAXSTORE);
        Chain->theta = doublevector(np);
        Chain->logP = doublevector(np);
        Chain->theta_store = (double**) malloc(np*sizeof(double*));
        memory += MAXSTORE*sizeof(int)+(2*MAXSTORE+np+np*np)*sizeof(double)+
        np*sizeof(double*);
      }

      /* --------------------------------------------- initialise parameter values */
      if (v == 2) printf("read parameters\n");
      Parameters(nt, nr, Ind->Chain[0].theta, Q);

      memory += Ind->chains*(sizeof(double)+sizeof(struct chain));

      /* ---------------------------------- storage for the posterior distribution */
      for (k = 0; k < np; k++) if (Use(FITTED_AUXILLARY, Q->P[k].type))
  Ind->Chain[0].theta_store[k] = doublevector(MAXSTORE);
      
      if (v == 2) printf("remove block\n");
      /* --------- remove any blocked parameters that are not declared using Par() */
      /* ---------------------------------------------- or are constant or derived */
      if (B->pblock != NULL)
  for (block = 0; block < B->nblocks; block++) {
    pb = B->pblock[block];
    npars = pb[0];
    for (i = 1; i <= pb[0]; i++) {
      if (pb[i] == -1) break;
      if (pb[i] >= Q->parameters       || 
    Q->P[pb[i]].type == NOTUSED || 
    Q->P[pb[i]].type == CONSTANT || 
    Q->P[pb[i]].type == INDIVIDUAL_CONSTANT || 
    Q->P[pb[i]].type == AUXILLARY) {
        fprintf(stderr, "Removing blocked parameter %d from treatment %d, "
          "individual %d as it is not declared\n", pb[i], nt, nr);
        for (j = i; j < npars; j++)
    pb[j] = pb[j+1];
        pb[npars] = -1;
        npars--;
        i--;
      }
    }
    if (npars == 0) {
      fprintf(stderr, "Block %d has no members\n", block);
      Error("");
    }
    pb[0] = npars;
  }

      if (v == 2) printf("add block\n");
      /* --------------------------- add any nonblocked parameters to extra blocks */
      nonblocked = integervector(np);
      for (k = 0; k < np; k++) if (Use(FITTED, Q->P[k].type)) 
  nonblocked[k] = TRUE;

      if (B->pblock != NULL)
  for (block = 0; block < B->nblocks; block++)
    for (i = 1; i <= B->pblock[block][0]; i++)
      nonblocked[B->pblock[block][i]] = FALSE;

      j = 0;
      for (k = 0; k < np; k++)
  if (nonblocked[k]) j++;

      if (j) {
  i = 0;
  newblock = TRUE;
  for (k = 0; k < np; k++) {
    if (newblock) {
      B->pblock = (int**) realloc(B->pblock, (B->nblocks+1)*sizeof(int*));
      B->pblock[B->nblocks] = integervector(Sim->parsperblock+1);
      newblock = FALSE;
    }
    // TODO CHANGE SO both can be blocked toegther
    if (nonblocked[k] && (Q->P[k].number != CATEGORY ||
        Q->P[k].hashyper == FALSE)) {
      nonblocked[k] = FALSE;
      j--;
      B->pblock[B->nblocks][i+1] = k;
      if (++i == Sim->parsperblock) {
        B->pblock[B->nblocks][0] = Sim->parsperblock;
        i = 0;
        newblock = TRUE;
        B->nblocks++;
      }
    }
  }
  if (i == 0) 
    B->nblocks--;
  else
    B->pblock[B->nblocks][0] = i;
  B->nblocks++;
  /* ----------------------------- block all categorical parameters together */
  if (j) {
    B->pblock = (int**) realloc(B->pblock, (B->nblocks+1)*sizeof(int*));
    B->pblock[B->nblocks] = integervector(j+1);
    i = 0;
    for (k = 0; k < np; k++)
      if (nonblocked[k] && Q->P[k].number == CATEGORY)
        B->pblock[B->nblocks][++i] = k;
    B->pblock[B->nblocks][0] = i;
    B->nblocks++;
  }    
  // need to check this works not sure need to do
  /* ---------------------------------- block parameters with hyper together */
  if (j) {
    B->pblock = (int**) realloc(B->pblock, (B->nblocks+1)*sizeof(int*));
    B->pblock[B->nblocks] = integervector(j+1);
    i = 0;
    for (k = 0; k < np; k++)
      if (nonblocked[k] && Q->P[k].hashyper)
        B->pblock[B->nblocks][++i] = k;
    B->pblock[B->nblocks][0] = i;
    B->nblocks++;
  }    
      }
      free(nonblocked);

      if (v == 2) printf("setup chain blocks\n");
      /* --------------------------------------------- setup blocks for each chain */
      for (c = 0; c < Ind->chains; c++) {
  Chain = &Ind->Chain[c];
  Chain->PB = (struct parameter_block*) 
    malloc(B->nblocks*sizeof(struct parameter_block));
  for (block = 0; block < B->nblocks; block++) {
    PB = &Chain->PB[block];
    PB->npar = B->pblock[block][0];
    PB->pblock = integervector(PB->npar);
    for (i = 0; i < PB->npar; i++)
      PB->pblock[i] = B->pblock[block][i+1];
    if (Q->P[PB->pblock[0]].distribution == CATEGORICAL)
      PB->categorical = TRUE;
    else {
      PB->categorical = FALSE;
      /* PB->L = doublematrix(PB->npar, PB->npar); */
      PB->L = gsl_matrix_alloc(PB->npar, PB->npar);
    }
  }
      }
    }
  //  memory += Data->Treatment[nt].nind*sizeof(struct sim_individual);

  /* ------------------------------------------- initialise hyper parameters store */
  if (HP->nhyper_parameters) {
    HP->hyper = doublevector(HP->nhyper_parameters);
    HP->hyper_store = doublematrix(HP->nhyper_parameters, MAXSTORE);
  }

  if (v == 2) printf("m3: %d\n", memory);

  memory += Sim->threads*sizeof(struct thread);
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      np = Sim->Treatment[nt].Individual[nr].Q.parameters;
      for (k = 0; k < np; k++) 
  if (Sim->Treatment[nt].Individual[nr].Q.P[k].type)
    memory += MAXSTORE*sizeof(double);
      }
  if (Sim->verbose) printf("memory = %.1fMB\n", memory/1048576.);
  
  /* -------------------------------------------------------------- output logfile */
  fp = fopen(Sim->logfile, "a");
  fprintf(fp, "random seed=%d\n", Sim->seed);
  fprintf(fp, "memory = %.1fMB\n", memory/1048576.);
  fprintf(fp, "individual parameters=%d\n", Sim->maxnp);
  fprintf(fp, "hyper parameters=%d\n", HP->nhyper_parameters);

  fprintf(fp, "Parameter blocks: Treatment Individual #parameters\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Chain = &Ind->Chain[0];
      for (block = 0; block < Ind->Q.Blocks.nblocks; block++) {
  PB = &Chain->PB[block];
  fprintf(fp, "%d: %d %d %d: ", block, nt, nr, PB->npar);
  for (i = 0; i < PB->npar; i++)
    fprintf(fp, " %d", PB->pblock[i]);
  fprintf(fp, "\n");
      }
    }

  fprintf(fp, "variable\tname\n");
  for (i = 0; i < Sim->V.variables; i++)
    fprintf(fp, "%d\t\"%s\"\n", i, Sim->V.W[i].name);

  fprintf(fp, "constant\tname\tnumber\tvalue\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Q = &Ind->Q;
      fprintf(fp, "replicate %d %d\n", nt, nr);
      for (i = 0; i < Q->parameters; i++) {
  if (Use(CONSTANT, Q->P[i].type))
    fprintf(fp, "%d\t\t\"%s\"\t%d\t%.3e\n", i,
      Q->P[i].name, Q->P[i].number, Q->P[i].location);
  else if (Use(INDIVIDUAL_CONSTANT, Q->P[i].type))
    fprintf(fp, "%d\t\t\"%s\"\t%d\tvariable\n", i, Q->P[i].name,
      Q->P[i].number);
      }
    }

  fprintf(fp, "parameter\thyperparameters\tname\tnumber\ttype\tdistribution\t"
    "location\tscale\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Q = &Ind->Q;
      fprintf(fp, "replicate %d %d\n", nt, nr);
      for (i = 0; i < Q->parameters; i++) 
  if (Use(FITTED_DERIVED_AUXILLARY, Q->P[i].type)) {
    fprintf(fp, "%d\t\t", i);
    if (Q->P[i].hashyper) fprintf(fp, "yes");
    fprintf(fp, "\t\t\"%s\"\t%d\t%d\t%d\t", Q->P[i].name, Q->P[i].number, 
      Q->P[i].type, Q->P[i].distribution);
    if (!Q->P[i].hashyper) 
      if (Q->P[i].distribution == BINOMIAL ||
    Q->P[i].distribution == GAMMA ||
    Q->P[i].distribution == INVGAMMA ||
    Q->P[i].distribution == BETA)
        fprintf(fp, "%.3e\t%.3e", Q->P[i].shape, Q->P[i].scale);
      else
        fprintf(fp, "%.3e\t%.3e", Q->P[i].location, Q->P[i].scale);
    else if (Q->P[i].distribution == BINOMIAL)
      fprintf(fp, "%.3e", Q->P[i].shape);
    fprintf(fp, "\n");
  }
    }

  fprintf(fp, "hyper parameters\n");
  fprintf(fp, "parameter\tname\tnumber\thyper distribution\tlocation\tscale\n");
  for (i = 0; i < HP->nhyper_parameters; i++) if (H[i].use)
    fprintf(fp, "%d\t\"%s\"\t%d\t%d\t%.3e\t%.3e\n", i, H[i].name, 
      H[i].number, H[i].hdist, H[i].location, H[i].scale);

  fprintf(fp, "dof\n");
  /*
  fprintf(fp, "hyperblocks\n");
  for (i = 0; i < HP->nblocks; i++) {
    HB = &HP->block[i];
    printf("%d %d\n", i, HB->use);
    if (HB->use) {
      fprintf(fp, "BLOCK %d\n", i);
      fprintf(fp, "nhypers = %d\n\t", HB->nhypers_used);
      for (j = 0; j < HB->nhypers_used; j++)
    fprintf(fp, "%d ", HB->hyper_used[j]);
      fprintf(fp, "\n");
      fprintf(fp, "associated pars = %d\n", HB->nassociated_pars);
      for (j = 0; j < HB->nassociated_pars; j++) {
    fprintf(fp, "\tpar %d: dependent block=%d ntreatments=%d\n\t",
      HB->associated_par[j], HB->dependent_block[j],
      HB->ntreatments_with_parameter[j]);
    for (k = 0; k < HB->ntreatments_with_parameter[j]; k++)
      fprintf(fp, "%d ", HB->treatment_with_parameter[j][k]);
    fprintf(fp, "\n");

    fprintf(fp, "\thypers with parameter, treatment combination\n");
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
        fprintf(fp, "\t\ttrt %d: ", nt);
        for (k = 0; k < HB->nhypers_with_treatment[j][nt]; k++)
          fprintf(fp, "%d ", HB->hyper_with_treatment[j][nt][k]);
        fprintf(fp, "\n");
    }
      }
    }
  }
  */
  fprintf(fp, "trt\tind\tpars\tparameters\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      fprintf(fp, "%d %d\t%d\t", nt, nr, Q->npar);
      for (k = 0; k < Q->parameters; k++) 
  if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) 
    fprintf(fp, "%d ", k);
      fprintf(fp, "\n");
    }

  fprintf(fp, "factors: %s\n", Sim->factors);

  fprintf(fp, "trt ind levels\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
      fprintf(fp, "%d %d %s\n", nt, nr, Data->Treatment[nt].levels);

  fprintf(fp, "posterior samples=%d\n", Sim->samples);
  fprintf(fp, "thinning=%d\n", Sim->thin);
  fclose(fp);

  if (!Sim->mcmc) {
    printf("mcmc = 0: not running Markov chain\n");
    return;
  }

  /* ------------------------------------------------- initialise parameter values */
  /* ------- need to do this before Saturated Likelihood for measurement variances */
  if (Sim->verbose) printf("Initialising parameter values\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      ReadParameters(Sim, &Ind->Q, nt, nr, Ind->Chain[0].theta);
    }

  /* --------------------------------------------- calculate saturated likelihoods */
  if (V->variables || finite(Sim->saturatedlogL))
    SaturatedLikelihood(Sim, Data, V->variables);

  /* --------------------------------------------- create a new temperature ladder */
  if (Sim->population) {
    v = Sim->verbose;
    if (v != 2) Sim->verbose = 0;
    mintemp = INFINITY;

    /* -------------------------------------- first pass: find minimum temperature */
    redotemp = (int**) malloc(Data->ntrt*sizeof(int*));
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
      redotemp[nt] = integervector(Data->Treatment[nt].nind);
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  if (Ind->chains <= 2) continue;
  Chain  = &Ind->Chain[0];

  /* -- temporarily set temperature to 0 to sample the prior to find mintemp */
  i = Ind->chains;
  Ind->chains = 1;
  Ind->temp[0] = 0;
  Chain->minlogL = INFINITY;
  RunMcMC(Sim, Data, nt, nr, MINTEMP);
  Ind->chains = i;
  Ind->temp[0] = 1;
  for (block = 0; block < Ind->Q.Blocks.nblocks; block++)
    redotemp[nt][nr] = redotemp[nt][nr] | Chain->PB[block].error;

  /* ---------------------------------- construct ladder if sampling worked */
  if (!redotemp[nt][nr]) {
    if (Chain->minlogL < mintemp) mintemp = Chain->minlogL;
    if (HP->nhyper_parameters == 0)
      ConstructLadder(Sim, Ind, fabs(1./Chain->minlogL), nt, nr);
  }
      }
    }

    /* ---- second pass using global minimum temperature for all individuals that */
    /* ------------------ didn't work first time or if there are hyper parameters */
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  if (redotemp[nt][nr] || HP->nhyper_parameters) 
    ConstructLadder(Sim, Ind, fabs(1./mintemp), nt, nr);
      }
    if (HP->nhyper_parameters) {
      nt = 0;
      while (!Sim->trtused[nt]) nt++;
      nr = 0;
      while (!Sim->indused[nt][nr]) nr++;
      Sim->chains = Sim->Treatment[nt].Individual[nr].chains;
    }
    Sim->verbose = v;
    for (nt = nt1; nt <= NT2; nt++)
      if (Sim->trtused[nt])
  free(redotemp[nt]);
    free(redotemp);
  }

  /* --------- check prior densities ------------*/
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      InitialiseIndividualCov(Sim, Data, Ind, nt, nr, ADAPTIVE_INDIVIDUAL);
    }

  /* ------------------- setup covariance matrices for individual level parameters */
  /* ------------ memory constraints means that each individual is done seperately */
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
      RunMcMC(Sim, Data, nt, nr, ADAPTIVE_INDIVIDUAL);

  /* ---------------------------------------------------------- non-adaptive phase */
  MCMC(Sim, Data, NONADAPTIVE, Sim->phase2adjust, ALLTRT, ALLIND);
  
  /* ------------------------- output hyperparameter posterior means and precision */
  /* if (Sim->hypers_used) { */
    /* double tmpm, tmpv; */
    /* sprintf(filename, "Results/%s.hypm", Sim->outexpno); */
    /* sprintf(backup, "cp %s %s.bak", filename, filename); */
    /* system(backup); */
    /* fp = fopen(filename, "w"); */
    /* HS = &Sim->H; */
    /* for (k = 0; k < HS->npar; k++) { */
    /*   meanvar(HS->hypermean_store[k], Sim->samples, &tmpm, &tmpv); */
    /*   fprintf(fp, "-1 %d %.3e\n", HS->par_indices[k], tmpm); */
    /* } */
    /* for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) { */
    /*   HT = &Sim->Treatment[nt].H; */
    /*   for (k = 0; k < HT->npar; k++) { */
    /*   meanvar(HT->hypermean_store[k], Sim->samples, &tmpm, &tmpv); */
    /*   fprintf(fp, "%d %d %.3e\n", nt, HT->par_indices[k], tmpm); */
    /*   meanvar(HT->hyperprec_store[k], Sim->samples, &tmpm, &tmpv); */
    /*   fprintf(fp, "%d %d %.3e\n", nt, HT->par_indices[k], tmpm); */
    /*   } */
    /* } */
    /* fclose(fp); */
  /* } */

  /* --------------------------------------- output hyperparameter posterior sample*/
  
  if (HP->nhyper_parameters) {
    sprintf(filename, "Results/%s.hypr", Sim->outexpno);
    sprintf(backup, "cp %s %s.bak", filename, filename);
    system(backup);
    fp = fopen(filename, "w");
    fprintf(fp, "samples %d\n", Sim->samples);
    for (k = 0; k < HP->nhyper_parameters; k++) if (H[k].use) { 
      fprintf(fp, "%d ", k);
      for (i = 0; i < Sim->samples; i++)
  fprintf(fp, "%.3e ", HP->hyper_store[k][i]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* ------------------------------------------- output parameter posterior sample */
  sprintf(filename, "Results/%s.post", Sim->expno);
  sprintf(backup, "cp %s %s.bak", filename, filename);
  system(backup);

  sprintf(filename,"Results/%s.post", Sim->outexpno);
  fp = fopen(filename, "w");
  if (Sim->verbose) printf("Discarding %d samples\n", (int)(Sim->discard/Sim->thin));
  fprintf(fp, "samples %d\n", Sim->samples);
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt]) {
    T = &Sim->Treatment[nt];
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      Dind = &Data->Treatment[nt].Individual[nr];
      Chain0 = &Ind->Chain[0];
      Q = &Ind->Q;

      fprintf(fp, "individual %d %d %d ", nt, nr, Q->npar);
      for (k = 0; k < Q->parameters; k++) 
  if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) 
    fprintf(fp, "%d ", k);
      fprintf(fp, "\n");
      for (i = 0; i < Sim->samples; i++) {
  for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED_AUXILLARY, Q->P[k].type)) {
      if (Q->P[k].number == REAL)
        tmptheta[k] = Chain0->theta_store[k][i];
      else
        tmptheta[k] = rint(Chain0->theta_store[k][i]);
    }
    else if (Use(CONSTANT, Q->P[k].type))
      tmptheta[k] = Q->P[k].location;
    else if (Use(INDIVIDUAL_CONSTANT, Q->P[k].type))
      tmptheta[k] = Q->P[k].location;

  TP.mode = OUTPUT;
  TP.n = Dind->np;
  TP.i = 0;
  TP.t = Dind->Y;
  // TP.t = doublevector(TP.n);
  // for (k = 0; k < TP.n; k++)
  //   TP.t[k] = Dind->Y[k][0];

  k = nsteps(Sim, Dind);
  y = doublematrix(k, 1+Sim->V.variables);
  Clear(y, k, Sim->V.variables);
  DerivedParameters(nt, nr, k, tmptheta, y, &TP);
  // free(TP.t); 
  free(y[0]); free(y);

  for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED_DERIVED_AUXILLARY, Q->P[k].type)) {
      if (Q->P[k].number == REAL)
        fprintf(fp, "%.4e ", tmptheta[k]);
      else
        fprintf(fp, "%ld ", lrint(tmptheta[k]));
    }
  fprintf(fp, "%.3e %.3e %d\n", Chain0->logpoststore[i], Chain0->logLstore[i],
    Ind->chainno[i]);
      }
    }
  }
  fclose(fp);

  /* ---------------------------------------------- output estimate of mode (MAP) */
  sprintf(filename, "Results/%s.best", Sim->outexpno);
  sprintf(backup, "cp %s %s.bak", filename, filename);
  system(backup);
  fp = fopen(filename, "w");

  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Q = &Sim->Treatment[nt].Individual[nr].Q;
      for (k = 0; k < Q->parameters; k++) 
  if (Use(FITTED, Q->P[k].type) || Use(INDIVIDUAL_CONSTANT, Q->P[k].type)) {
    fprintf(fp, "%d %d %d ", nt, nr, k);
    if (Q->P[k].number == REAL)
      fprintf(fp, "%e\n", Sim->Treatment[nt].Individual[nr].Chain[0].theta[k]);
    else
      fprintf(fp, "%ld\n", 
        lrint(Sim->Treatment[nt].Individual[nr].Chain[0].theta[k]));
  }
    }
  fclose(fp);

  
  /* ------------------------------------------------- output marginal likelihoods */
  sprintf(filename, "Results/%s.diag", Sim->outexpno);
  fp = fopen(filename, "w");
  fprintf(fp, "Marginal log-Likelihoods\n");
  if (Sim->population) {
    if (Sim->xrate) {
      sprintf(backup, "cp %s %s.bak", Sim->tempfile, Sim->tempfile);
      system(backup);
      fpl = fopen(Sim->tempfile, "w");
    }
    for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  fprintf(fp, "%d %d ", nt, nr);
  mean = doublevector(Ind->chains);
  var = doublevector(Ind->chains);
  hmean = doublevector(Ind->chains);
  error2 = doublevector(Ind->chains);

  if (Sim->xrate) {
    ElogL.n = Ind->chains;
    ElogL.x = doublevector(ElogL.n);
    ElogL.y = doublevector(ElogL.n);
    ElogL.m = doublevector(ElogL.n);
    ElogL.c = doublevector(ElogL.n);
  }
  for (c = 0; c < Ind->chains; c++) {
    lnL = Ind->Chain[c].logLstore;
    nadj = r = Sim->samples;

    /* ----------------------- order log-likelihoods from highest to lowest */
    HeapSort(lnL-1, r);
    for (i = 0; i < r/2-1; i++) {
      f = lnL[i];
      lnL[i] = lnL[r-1-i];
      lnL[r-1-i] = f;
    }

    /* ----------------- find a measure of the scale of the log-likelihoods */
    f = lnL[0]-lnL[r/2];

    /* ----------------------------- remove extremely small log-likelihoods */
    if (f != 0) while(lnL[nadj-1] < lnL[r/2]-1000*f) nadj--;
    
    /* -------------------------------------- calculate mean log-likelihood */
    meanvar(lnL, nadj, &mean[c], &var[c]);
    error2[c] = var[c]/(double)nadj;

    if (Sim->xrate) {
      ElogL.x[c] = log(Ind->temp[c]);
      ElogL.y[c] = mean[c];
    }
  }
  /* ---------------------------------------------- log-marginal likelihood */
  logML = 0;
  for (c = 0; c < Ind->chains-1; c++)
    logML += (Ind->temp[c]-Ind->temp[c+1])*(mean[c]+mean[c+1])/2.;

  /* --------------------------- error, small compared to monte carlo error */
  seML = 0;
  for (c = 1; c <= Ind->chains-3; c++)
          seML += pow(Ind->temp[c+1]-Ind->temp[c-1], 2)*error2[c];
        seML += pow(Ind->temp[1]-Ind->temp[0], 2)*error2[0];
        seML += pow(Ind->temp[Ind->chains-2]-Ind->temp[Ind->chains-3], 2)*
          error2[Ind->chains-2];
        seML /= 2.*(double)nadj;

  fprintf(fp, "%e %e\n", logML, sqrt(seML));

  if (Sim->xrate) {
    c = 0;
    ConstructPiecewiseLine(&ElogL);
    tmpt[c] = 1;
    do {
      c++;
      ReturnPiecewiseLine(log(tmpt[c-1]), &f, &df, &ElogL);
      tmpt[c] = tmpt[c-1] - sqrt(-log(Sim->xrate)*tmpt[c-1]/df);
    } while (tmpt[c] > Ind->temp[Ind->chains-2]);
    fprintf(fpl, "%d %d %d\n", nt, nr, c);
    for (i = 0; i < c; i++)
      fprintf(fpl, "%e\n", tmpt[i]);

    free(ElogL.x); free(ElogL.y); free(ElogL.m); free(ElogL.c);
  }
  free(var); free(mean); free(hmean);
      }
    if (Sim->xrate) fclose(fpl);
  }
    
  /* ------------------------------------------------------------- output deviance */
  fprintf(fp, "Deviance [max(log(posterior)), max(log(likelihood)), "
    "posterior mean deviance, deviance, Ndata]\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      double tmpm, tmpv;
      Ind = &Sim->Treatment[nt].Individual[nr];
      meanvar(Ind->Chain[0].logLstore, Sim->samples, &tmpm, &tmpv);

      Chain0 = &Sim->Treatment[nt].Individual[nr].Chain[0];
      Dind = &Data->Treatment[nt].Individual[nr];
      fprintf(fp, "%d %d %e %e %e %e %d\n", nt, nr, Chain0->logpost, 
        Chain0->logmaxL, tmpm, -2.*(Chain0->logmaxL-Dind->saturated),
        Dind->n);
    }

  /* --------------------------------------------- output Gelman Rubin convergence */
  fprintf(fp, "Gelman Rubin statistics\n");
  GelmanRubinStats(Sim, Data, fp);

  /* ----------------------------------------------- output density exchange rates */
  fprintf(fp, "Density exchange rates\n");
  minX = 1;
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      fprintf(fp, "%d %d\n", nt, nr);
      for (c = 0; c < Ind->chains-1; c++) {
  X = Ind->Chain[c].exchangeaccept/Ind->Chain[c].exchangeproposed;
  fprintf(fp, "%d %e %f\n", c, Ind->temp[c], X);
  if (X < minX) minX = X;
      }
    }
  if (minX < 0.5) 
    fprintf(stderr, "WARNING! Minimum between-density exchange rate is "
      "less than 50%%. Use 'bag x' to examine exchange rates.\n");
  fclose(fp);

  /* --------------------------------------------------- output chain temperatures */
  sprintf(filename, "Results/%s.temp", Sim->outexpno);
  sprintf(backup, "cp %s %s.bak", filename, filename);
  system(backup);
  fp = fopen(filename, "w");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      fprintf(fp, "%d %d %d\n", nt, nr, Ind->chains-1);
      for (c = 0; c < Ind->chains-1; c++)
  fprintf(fp, "%e\n", Ind->temp[c]);
    }
  fclose(fp);
}

void ConstructLadder(Simulation *Sim, struct sim_individual *Ind, 
         double mintemp, int nt, int nr)
{
  int c = 0;

  if (mintemp < Ind->temp[Ind->chains-2])
    fprintf(stderr, "treatment %d, individual %d: temperature %e is below the "
      "lowest ladder temperature of %e\n", nt, nr, mintemp,
      Ind->temp[Ind->chains-2]);

  while (c < Ind->chains-1 && Ind->temp[c++] > mintemp);
  Ind->temp[0]   = 1;
  Ind->temp[c-1] = mintemp;
  Ind->temp[c]   = 0;
  Ind->chains    = c+1;
}

void RunMcMC(Simulation* Sim, struct data* Data, int nt, int nr, int phase)
{
  int block, c;
  struct parameter_block *PB;
  struct sim_individual *Ind = &Sim->Treatment[nt].Individual[nr];

  for (c = 0; c < Ind->chains; c++)
    for (block = 0; block < Ind->Q.Blocks.nblocks; block++) {
      PB = &Ind->Chain[c].PB[block];
      PB->theta_store = doublematrix(PB->npar, MAXSTORE);
    }

  if (InitialiseIndividualCov(Sim, Data, Ind, nt, nr, phase)) Error("");
  MCMC(Sim, Data, phase, Sim->phase1adjust, nt, nr);

  for (c = 0; c < Ind->chains; c++)
    for (block = 0; block < Ind->Q.Blocks.nblocks; block++) {
      PB = &Ind->Chain[c].PB[block];
      free(PB->theta_store[0]);
      free(PB->theta_store);
    }
}

int Check(Simulation *Sim, int k, double logP, double theta, int nt, int nr, 
    struct parameter *P, double scale) 
{
  /* -------------------------------- check parameter has a valid distribution */
  if (logP == M_PI) {
    fprintf(stderr, "parameter %d has no prior distribution\n", k);
    Error("");
  }
  /* -------------------------------- check parameter is within range of prior */
  if (isinf(logP) || isnan(logP)) {
    fprintf(stderr, "parameter %d=%e is out of prior range for "
      "treatment %d, individual %d (logP=%f)\n", k, theta, nt, nr, logP);
    Error("");
  }
  /*
    y=ln(P)+0.5*ln(s^2)
    y=ln(P * s)

    P=(2*pi*s^2)^(-1/2)*exp(-0.5*((x-mu)/s)^2)
    ln(P)=-0.5*ln(2*pi*s^2)-0.5*((x-mu)/s)^2
    x-mu=3s
    ln(P)=-0.5*ln(2*pi*s^2)--0.5*(3s/s)^2
    ln(P)=-0.5*ln(2*pi*s^2)--0.5*(3)^2
    y = -0.5*ln(2*pi*s^2)-0.5*(3)^2 + 0.5*ln(s^2)
    NORMAL: -0.5*(3)^2 = -4.5
    
   */
  
  /* ---------------------------- check parameter value is close to prior mean */
  /* if (!P->hashyper && P->distribution != CATEGORICAL &&  */
  /*     P->distribution != USERDEFINED && P->distribution != UNIFORM && */
  /*     logP + 0.5*log(scale) < Sim->minprior) { */
  /*   fprintf(stderr, "prior density of parameter %d=%e is very low " */
  /*       "for treatment %d, individual %d\n", k, theta, nt, nr); */
  /*   fprintf(stderr, "logP=%f scale=%f\n", logP, scale); */
  /*   if (!Sim->ignorewarning) Error(""); */
  /* } */
  /* else if (P->hashyper && (l=logP - 0.5*log(scale) < Sim->minprior)) { */
  /*   printf("%d %d %d %f %f %f\n", k, nt, nr, logP, theta, scale);  */
  /*   if (!Sim->ignorewarning) Error("hyper scale"); */
    //    Error("hyper scale");
    /* i = H->par_indices[k]; */
    /* fprintf(stderr, "prior density of hyper parameter %d=%e is very low " */
    /*       "for treatment %d, individual %d\n", i, theta, nt, nr); */
    /* fprintf(stderr, "l=logP + 0.5*log(variance) < minprior\n"); */
    /* fprintf(stderr, "l=%e\nminprior=%e\nlogP=%e\n0.5*log(variance)=%e\n", */
    /*       l, Sim->minprior, logP, -0.5*log(H->hyperprec.a[k+H->npar*k])); */
    /* fprintf(stderr, "perhaps set usehypermeans 1 in pars files\n"); */
    /* fprintf(stderr, "p[%d] = %e;\n", i, H->hypermean.a[k]); */
    /* Error(""); */
  /* } */
  return TRUE;
}

int InitialiseIndividualCov(Simulation* Sim, struct data* Data, 
          struct sim_individual *Ind, int nt, int nr, int phase)
{
  int i, k, c, block, flag, status = SUCCESS;
  struct chain *Chain, *Chain0;
  struct parameter_block *PB, *PB0;
  struct parameterset *Q = &Ind->Q;
  struct parameter *P;
  Hyperblock *HB = Sim->HP.block;
  double *hyper = Sim->HP.hyper;
  double *logP, *theta, sumlogP, l;

  Chain0 = &Ind->Chain[0];
  logP = Chain0->logP;
  theta = Chain0->theta;

  /* -------------------------------------------------------------- get parameters */
  ReadParameters(Sim, Q, nt, nr, theta);

  /* ----------------------------------------------------------------- error check */
  for (k = 0; k < Q->parameters; k++) if (Use(FITTED, Q->P[k].type))
    logP[k] = M_PI;
  sumlogP = sumlogPriors(Sim, nt, theta, logP, Q, hyper);

  for (k = 0; k < Q->parameters; k++) {
    P = &Q->P[k];
    if (Use(FITTED, P->type)) {
      if (P->hashyper)
  Check(Sim, k, logP[k], theta[k], nt, nr, P,
        Hypervalue(&HB[P->hyper_scale_block_index], hyper, P->index, nt));
      else
  Check(Sim, k, logP[k], theta[k], nt, nr, P, P->scale);
    }
  }

  /* ------------------------------------------- check likelihood is a real number */

  Chain0->logL = FulllogLikelihood(Sim, Data, Chain0->theta, nt, nr, &status);
  Chain0->logpost = Ind->temp[0]*Chain0->logL + sumlogP;
  Chain0->logmaxL = Chain0->logL;

  flag = FALSE;
  if (isinf(Chain0->logL)) {
    fprintf(stderr, "likelihood is infinite for treatment %d, individual %d\n",
            nt, nr);
    if (status) fprintf(stderr, "Model returned error %d\n", status);
    flag = TRUE;
  }
  if (flag) Error("");
  if (isnan(Chain0->logL)) {
    fprintf(stderr, "likelihood is NaN for treatment %d, individual %d\n"
      "you may be taking the log of a zero or a negative number", nt, nr);
    flag = TRUE;
  }
  if (flag) Error("");

  /* ---------------------------------------------- initialise covariance matrices */
  for (block = 0; block < Ind->Q.Blocks.nblocks; block++) {
    PB = &Chain0->PB[block];
    if (!PB->categorical) {
      PB->inrange  = 0;
      PB->covscale = 1;
      PB->covstart = 0;
      PB->error = FALSE;
      PB->covfinal = FALSE;

      gsl_matrix_set_zero(PB->L);

      for (i = 0; i < PB->npar; i++) {
  k = PB->pblock[i];
  l = Q->P[k].number == REAL ? 
    theta[k]/Sim->initscale : MAX(1, theta[k]/Sim->initscale);
  gsl_matrix_set(PB->L, i, i, l);
      }
    }
  }

  /* -------------------------------------------------- copy chain 0 to all chains */
  for (c = 1; c < Ind->chains; c++) {
    Chain = &Ind->Chain[c];

    Chain->logL    = Chain0->logL;
    Chain->logpost = Ind->temp[c]*Chain->logL + sumlogP;
    for (k = 0; k < Q->parameters; k++) {
      if (Use(FITTED_CONSTANT, Q->P[k].type))
  Chain->theta[k] = theta[k];
      if (Use(FITTED, Q->P[k].type))
  Chain->logP[k]  = logP[k];
    }

    for (block = 0; block < Ind->Q.Blocks.nblocks; block++) {
      PB = &Chain->PB[block];
      if (!PB->categorical) {
  PB0 = &Chain0->PB[block];
  PB->inrange  = PB0->inrange;
  PB->covscale = PB0->covscale;
  PB->covstart = PB0->covstart;
  PB->error    = PB0->error;
  PB->covfinal = PB0->covfinal;

  gsl_matrix_memcpy(PB->L, PB0->L);
      }
    }
  }
  return flag;
}

void MCMC(Simulation* Sim, struct data* Data, int phase, int adjust, 
    int ntp, int nrp)
{
  int i, b, c, cs, thread, k, nt, nr, change, mindensity, anychange;
  int v = Sim->verbose, it = -1, *globaln = &it, n = -1;
  double rate, tmp, DlogP, dt, prop_logL;
  double exchangerate, *store_theta, *store_int_theta, *store_logP;
  double minexchangerate, *hyper;
  struct sim_treatment *T;
  struct sim_individual *Ind;
  struct thread Thread;
  struct chain *Chain, *ChainT, *Chain0, *Chain1;
  struct data_individual *D;
  struct parameterset *Q;
  struct parameter_block *PB;
  gsl_rng *stream;
  FILE *fp;
  static int firsttime = TRUE;
  
  if (firsttime) {
    fp = fopen("/tmp/bayes.log", "w");
    firsttime = FALSE;
  }
  else
    fp = fopen("/tmp/bayes.log", "a");

  /* --------------------------------------------------------- reset chain storage */
  Sim->storeno = 0;
  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
    T = &Sim->Treatment[nt];
    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      for (c = 0; c < Ind->chains; c++) {
  Chain = &Ind->Chain[c];
  Chain->exchangeaccept = Chain->exchangeproposed = 0;
  Chain->accept = Chain->proposed = 0;
  Chain->done = FALSE;
  Chain->loop = FALSE;
  Chain->laststoreno = Chain->storeno;
  Chain->storeno = 0;
  for (b = 0; b < Ind->Q.Blocks.nblocks; b++) {
    PB = &Chain->PB[b];
    PB->accept = PB->proposed = 0;
    PB->storeno = 0;
    PB->covstart = 0;
    PB->covfinal = FALSE;
    PB->done = FALSE;
  }
      }
    }
  }
  if (v && phase == ADAPTIVE_INDIVIDUAL) {
    printf("-Adaptive-\n");
    fprintf(fp, "-Adaptive-\n");
  }

# pragma omp parallel for            \
  private(i, b, c, k, nt, nr, Thread, Chain, ChainT, Chain0, Chain1,  \
    T, Ind, change, rate, DlogP, dt, prop_logL, anychange,  \
    store_theta, store_logP, exchangerate, minexchangerate,  \
    cs, mindensity, stream, tmp, Q, hyper, PB,      \
    store_int_theta)            \
  firstprivate(n, phase, v, fp)
  
  for (thread = 0; thread < Sim->threads; thread++) {
    SetupThread(&Thread, Sim, Data, phase, ntp, nrp, thread);

    stream = Thread.stream;
    hyper = Thread.hyper;
    store_int_theta = doublevector(Sim->maxnp);
    store_theta = doublevector(Sim->maxnp);
    store_logP = doublevector(Sim->maxnp);

    while (TRUE) {
      n++;
      if (Done(Sim, Data, ntp, nrp, phase)) {
  ClearThread(&Thread, Sim, Data, ntp, nrp);
  break;
      }

      /* ------------------------------------------ output chain descriptive stats */
#     pragma omp critical (covariance)
      {
  (*globaln)++;
  if (*globaln % adjust == 0) {
    if ((v && phase == NONADAPTIVE) || v == 2) {
      printf("\n+++++++++++++++ %.2f\n", GelmanRubinStats(Sim, Data, NULL));
      fprintf(fp, "\n+++++++++++++++ %.2f\n", GelmanRubinStats(Sim, Data,
                     NULL));
    }
    for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
      if ((v && phase == NONADAPTIVE) || v == 2) {
        if (phase == MINTEMP) {
    printf("-MinTemp-- %d\n", *globaln);
    fprintf(fp, "-MinTemp-- %d\n", *globaln);
        }
        else if (phase == ADAPTIVE_INDIVIDUAL) {
    printf("-Adaptive- %d\n", *globaln);
    fprintf(fp, "-Adaptive- %d\n", *globaln);
        }
        else
    if (*globaln-Sim->discard < 0) {
      printf("-Burnin--- %d\n", *globaln-Sim->discard);
      fprintf(fp, "-Burnin--- %d\n", *globaln-Sim->discard);
    }
    else {
      printf("---------- %d\n", *globaln-Sim->discard);
      fprintf(fp, "---------- %d\n", *globaln-Sim->discard);
    }
      }
      /* else if (v && phase == ADAPTIVE_INDIVIDUAL) */
      /*   printf("-Adaptive- %d\n", *globaln); */

       for (nr = NRS; nr <= NRE; nr++)  if (Sim->indused[nt][nr]) {
        Ind = &Sim->Treatment[nt].Individual[nr];
        D = &Data->Treatment[nt].Individual[nr];

        minexchangerate = 1;
        mindensity = 0;
        for (c = 0; c < Ind->chains; c++) {
    Chain = &Ind->Chain[c];
    if (Chain->exchangeproposed) {
      exchangerate = Chain->exchangeaccept/Chain->exchangeproposed;
      if (exchangerate < minexchangerate) {
        minexchangerate = exchangerate;
        mindensity = c;
      }
    }
        }
        if (phase == NONADAPTIVE && minexchangerate < Sim->minexchangerate) 
    v = HIGH;

        if (v) {
    for (c = 0; c < (v == 1 ? 1 : Ind->chains); c++) {
      if (phase == MINTEMP || phase == ADAPTIVE_INDIVIDUAL) {
        Chain = &Ind->Chain[c];
        for (b = 0; b < Ind->Q.Blocks.nblocks; b++) {
          PB = &Chain->PB[b];
          if (PB->done && c && b) continue;

          rate = PB->proposed ? PB->accept/PB->proposed : 0;
          if (v == HIGH)
      if (PB->categorical) {
        printf("%- 2u %- 2u | %- 2u %- 2u | categorical     ", 
         nt, nr, c, b);
        fprintf(fp, "%- 2u %- 2u | %- 2u %- 2u | categorical     ", 
          nt, nr, c, b);
      }
      else {
        printf("%- 2u %- 2u | %- 2u %- 2u |"
         " % 2lu % 4u % 4u %u %.1e ", 
         nt, nr, c, b, lrint(100.*rate), 
         PB->storeno, PB->covstart, PB->covfinal, 
         PB->covscale);
        fprintf(fp, "%- 2u %- 2u | %- 2u %- 2u |"
          " % 2lu % 4u % 4u %u %.1e ", 
          nt, nr, c, b, lrint(100.*rate), 
          PB->storeno, PB->covstart, PB->covfinal, 
          PB->covscale);
      }
          else if (v == LOW) {
      if (PB->categorical) {
        printf("%- 2u %- 2u | %- 2u | categorical     ", 
         nt, nr, b);
        fprintf(fp, "%- 2u %- 2u | %- 2u | categorical     ", 
         nt, nr, b);
      }
      else {
        printf("%- 2u %- 2u | %- 2u |"
         " % 2lu % 4u % 4u %u %.1e ", 
         nt, nr, b, lrint(100.*rate), 
         PB->storeno, PB->covstart, PB->covfinal, 
         PB->covscale);
        fprintf(fp, "%- 2u %- 2u | %- 2u |"
           " % 2lu % 4u % 4u %u %.1e ",
           nt, nr, b, lrint(100.*rate),
           PB->storeno, PB->covstart, PB->covfinal,
           PB->covscale);
      }
          }
          printf("| % 7.1f % 7.1f",
           Chain->logL, -2.*(Chain->logL - D->saturated)); 
          printf("\n"); 
          fprintf(fp, "| % 7.1f % 7.1f",
            Chain->logL, -2.*(Chain->logL - D->saturated)); 
          fprintf(fp, "\n"); 
        }
      }
      else if (phase == NONADAPTIVE) { 
        Chain = &Ind->Chain[c];
        exchangerate = Chain->exchangeproposed ?
          Chain->exchangeaccept/Chain->exchangeproposed : 0;
        if (v == HIGH) {
          printf("%- 2u %- 2u %- 2u | ", nt, nr, c);
          fprintf(fp, "%- 2u %- 2u %- 2u | ", nt, nr, c);
        }
        else {
          printf("%- 2u %- 2u | ", nt, nr);
          fprintf(fp, "%- 2u %- 2u | ", nt, nr);
        }
        for (b = 0; b < Ind->Q.Blocks.nblocks; b++) {
          PB = &Chain->PB[b];
          rate = PB->proposed ? PB->accept/PB->proposed : 0;
          if (b < Ind->Q.Blocks.nblocks-1) {
      printf("% 2lu ", lrint(100.*rate));
      fprintf(fp, "% 2lu ", lrint(100.*rate));
          }
          else {
      printf("% 2lu ", lrint(100.*rate));
      fprintf(fp, "% 2lu ", lrint(100.*rate));
          }
        }
        if (v == HIGH) {
          printf("| % 2lu %.2e %d | ", lrint(100.*exchangerate), 
           Ind->temp[c], Ind->chains);
          fprintf(fp, "| % 2lu %.2e %d | ", lrint(100.*exchangerate), 
           Ind->temp[c], Ind->chains);
        }
        else if (v == LOW && Ind->chains > 1) {
      printf("| % 2lu ", lrint(100.*minexchangerate));
      fprintf(fp, "| % 2lu ", lrint(100.*minexchangerate));
        }
        /*
          if (c == 0) { 
          for (k = 0; k < Sim->threads; k++) 
          tmplogL[k] = Sim->Thread[k].Treatment[nt]. \ 
          Individual[nr].Chain[c].logL;
          meanvar(tmplogL, Sim->threads, &mean, &var);
          printf("% 7.1f % 3.1f | % 7.1f % 7.1f", mean, sqrt(var),
          Chain->logL, -2.*(Chain->logL - D->saturated)); 
          }
        */
        printf("| % 7.1f % 7.1f % 7.1f", 
           Chain->logmaxL, Chain->logpost, 
           -2.*(Chain->logmaxL - D->saturated)); 
        printf("\n"); 
        fprintf(fp, "| % 7.1f % 7.1f % 7.1f",
           Chain->logL, Chain->logpost, 
           -2.*(Chain->logmaxL - D->saturated)); 
        fprintf(fp, "\n"); 
      }
    }
        }
        /* for (c = 0; c < Ind->chains; c++) { */
        /*   Chain = &Ind->Chain[c]; */
        /*   printf("%d %d %d %d\n", nt, nr, c, Chain->storeno); */
        /*   fprintf(fp, "%d %d %d %d\n", nt, nr, c, Chain->storeno); */
        /* } */
        if (phase == NONADAPTIVE && minexchangerate < Sim->minexchangerate) {
    fprintf(stderr, "density %d exchange rate = %.1f for %d,%d\n", 
      mindensity, 100.*minexchangerate, nt, nr);
    fprintf(fp, "density %d exchange rate = %.1f for %d,%d\n", 
      mindensity, 100.*minexchangerate, nt, nr);
    Error("Exchange rate too low: Try setting number of densities to "
          "1, or set minexchangerate to 0 in pars");
        }

        /* ----------------------------------- calculate covariance matrices */
        if (*globaln) {
    Q = &Sim->Treatment[nt].Individual[nr].Q;
    for (c = 0; c < Ind->chains; c++) {
      Chain = &Ind->Chain[c];
      for (b = 0; b < Q->Blocks.nblocks; b++) {
        PB = &Chain->PB[b];
        if ((phase == MINTEMP || phase == ADAPTIVE_INDIVIDUAL) && 
      !PB->categorical && !PB->done) 
          AdjustCovariance(Sim, PB, Q, phase);
        PB->accept = PB->proposed = 0;
      }
      Chain->accept = Chain->proposed = 0;
    }
        }
      }
    }
  }
  fflush(fp);
      }

      /* ------------------------------------------------- update hyper parameters */
      if (phase == NONADAPTIVE && Sim->HP.nhyper_parameters) {
  Gibbs(Sim, Data, &Thread, stream);

  /* ---------------------------------------- recalculate priors and logpost */
  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt])
    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Q = &Ind->Q;
      for (c = 0; c < Sim->chains; c++) {
        ChainT = &Thread.Treatment[nt].Individual[nr].Chain[c];
        ChainT->logpost = Ind->temp[c]*ChainT->logL + 
    sumlogPriors(Sim, nt, ChainT->theta, ChainT->logP, Q, hyper);
      }
    }

#       pragma omp critical (hyperstore)
  if (*globaln-Sim->discard > 0 && n%Sim->thin == 0) {
    for (k = 0; k < Sim->HP.nhyper_parameters; k++) if (Sim->HP.H[k].use)
        Sim->HP.hyper_store[k][Sim->storeno] = hyper[k];
    if (++Sim->storeno >= MAXSTORE) Sim->storeno = 0;
  }
      }

      /* -------------------------------------------- update individual parameters */
      for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
  for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
    Ind = &Sim->Treatment[nt].Individual[nr];
    Q = &Ind->Q;
    for (c = 0; c < Ind->chains; c++) {
      Chain  = &Ind->Chain[c];
      ChainT = &Thread.Treatment[nt].Individual[nr].Chain[c];

      /* ------------------------------------------- block update parameters */
      anychange = FALSE;
      for (b = 0; b < Q->Blocks.nblocks; b++) {
        PB = &Chain->PB[b];
#             pragma omp critical (proposed)
        PB->proposed++;

        /* ---------- save blocked parameters to be changed and their priors */
         for (i = 0; i < PB->npar; i++) {
    k = PB->pblock[i];
    store_theta[k] = ChainT->theta[k];
    store_logP[k]  = ChainT->logP[k];
        }

        /* --------------------------------------- save auxillary parameters */
        for (i = 0; i < Q->nauxillary; i++) {
    k = Q->auxillary[i];
    store_theta[k] = ChainT->theta[k];
        }
        
        Proposal(stream, ChainT->theta, PB);
    
        /* ----- store integer proposed parameters and convert with rint */
        for (i = 0; i < Q->ninteger; i++) {
    k = Q->integer[i];
    store_int_theta[k] = ChainT->theta[k];
    ChainT->theta[k] = rint(ChainT->theta[k]);
        }
    
        // calculate change in logP
        DlogP = 0;
        for (i = 0; i < PB->npar && finite(DlogP); i++) {
    k = PB->pblock[i];
    ChainT->logP[k] = logPrior(Sim, nt, ChainT->theta[k], &Q->P[k], 
             hyper);
    DlogP += ChainT->logP[k] - store_logP[k];
        }

        change = FALSE;
        if (finite(DlogP))
    Metropolis(Sim, nt, nr, ChainT, Ind->temp[c], DlogP,
         PB, phase, Data, stream, &prop_logL, &change);

        if (!change) {
    /* ------------------ restore blocked parameters if not accepted */
    for (i = 0; i < PB->npar; i++) {
      k = PB->pblock[i];
      ChainT->theta[k] = store_theta[k];
      ChainT->logP[k] = store_logP[k];
    }
    /* ---------------- restore auxillary parameters if not accepted */
    for (i = 0; i < Q->nauxillary; i++) {
      k = Q->auxillary[i];
      ChainT->theta[k] = store_theta[k];
    }
        } else {
    /*  restore real values of integer parameters values if accepted */
    for (i = 0; i < Q->ninteger; i++) {
      k = Q->integer[i];
      ChainT->theta[k] = store_int_theta[k];
    }
        }
        /* --------------------------- store samples for covariance matrix */
#             pragma omp critical (covstore)
        {
    if ((phase == MINTEMP || phase == ADAPTIVE_INDIVIDUAL) && 
        !PB->categorical && !PB->done && change) {
      for (k = 0; k < PB->npar; k++)
        PB->theta_store[k][PB->storeno] = ChainT->theta[PB->pblock[k]];
      
      if (phase == MINTEMP && prop_logL < Chain->minlogL) 
        Chain->minlogL = prop_logL;
      
      if (PB->covfinal) {
        PB->storeno++;
        if (PB->storeno == PB->covstart+Sim->covsamples) PB->done = TRUE;
      } else {
        PB->covstart++;
        PB->storeno++;
        if (PB->covstart > Sim->covsamples) {
          PB->covstart = Sim->covsamples;
          if (PB->storeno == PB->covstart) PB->storeno = 0;
        }
      }
    }
        }
        anychange = anychange | change;
      }

      /* -------------------------- save estimate of mode and max likelihood */
#           pragma omp critical (best)
      {
        if (ChainT->logpost > Chain->logpost)
    Best(Sim, Data, Q, ChainT, Chain, ntp, nrp, c, phase);
        if (c == 0 && ChainT->logL > Chain->logmaxL)
    Chain->logmaxL = ChainT->logL;
      }

      /* --------------------------------------- store samples for inference */
#           pragma omp critical (store) 
      {
        if ((phase == NONADAPTIVE && *globaln-Sim->discard > 0 && 
       n%Sim->thin == 0) ||
      (phase == ADAPTIVE_INDIVIDUAL && anychange)) {
    Chain->logLstore[Chain->storeno] = ChainT->logL;
    Chain->logpoststore[Chain->storeno] = ChainT->logpost;
    Ind->chainno[Chain->storeno] = omp_get_thread_num();
    
    if (c == 0) {
      for (k = 0; k < Q->parameters; k++) 
        if (Use(FITTED_AUXILLARY, Q->P[k].type))
          Chain->theta_store[k][Chain->storeno] = ChainT->theta[k];
        }
    if (++Chain->storeno >= MAXSTORE) {
      if (phase == ADAPTIVE_INDIVIDUAL) {
        Chain->storeno = 0;
        Chain->loop = TRUE;
      }
      else
        Error("here");
    }
    
    /* ------------------------------------------------- check if done */
    if (Chain->storeno == Sim->samples) Chain->done = TRUE;
        }
      }
    }
  }
      }

      /* --------------------------------------- exchange states between densities */
      if ((phase == ADAPTIVE_INDIVIDUAL || phase == NONADAPTIVE) && 
    Sim->population) {
  cs = uniformRNG(0, 1, stream) < 0.5? 0: 1;
  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt])
    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      for (c = cs; c < Ind->chains-1; c += 2) {
        Chain  = &Ind->Chain[c];
        Chain0 = &Thread.Treatment[nt].Individual[nr].Chain[c];
        Chain1 = &Thread.Treatment[nt].Individual[nr].Chain[c+1];
        dt = Ind->temp[c]-Ind->temp[c+1];
#             pragma omp critical (exchangeproposed)
        Chain->exchangeproposed++;
        
        if (Chain0->logL != Chain1->logL && 
      uniformRNG(0, 1, stream) < exp(dt*(Chain1->logL-Chain0->logL))) {
#               pragma omp critical (exchangeaccept)
    Chain->exchangeaccept++;

    tmp = Chain0->logpost;
    Chain0->logpost = Chain1->logpost + dt*Chain1->logL;
    Chain1->logpost = tmp - dt*Chain0->logL;
    
    tmp = Chain0->logL;
    Chain0->logL = Chain1->logL;
    Chain1->logL = tmp;
    
    for (k = 0; k < Ind->Q.parameters; k++) {
      if (Use(FITTED_AUXILLARY, Ind->Q.P[k].type)) {
        tmp = Chain0->theta[k];
        Chain0->theta[k] = Chain1->theta[k];
        Chain1->theta[k] = tmp;
      }
      if (Use(FITTED, Ind->Q.P[k].type)) {
        tmp = Chain0->logP[k];
        Chain0->logP[k] = Chain1->logP[k];
        Chain1->logP[k] = tmp;
      }
    }
        }
      }
    }
      }
    }
    free(store_int_theta);
    free(store_theta);
    free(store_logP);
  }
  fclose(fp);
}

void Metropolis(Simulation *Sim, int nt, int nr,
    struct chain *ChainT, double temp, double DlogP,
    struct parameter_block *PB, int phase, struct data *Data, 
    gsl_rng *stream, double *prop_logL, int *change)
{
  int status;
  double alpha;

  *prop_logL = FulllogLikelihood(Sim, Data, ChainT->theta, nt, nr, &status);
  alpha = temp*(*prop_logL-ChainT->logL) + DlogP;

  /* ---------------------------------------------------------- do Metropolis step */
  if (uniformRNG(0, 1, stream) < exp(alpha)) {
    ChainT->logpost += alpha;
    ChainT->logL = *prop_logL;
    *change = TRUE;
#   pragma omp critical (accept) 
    PB->accept++;
  }
}

double FulllogLikelihood(Simulation *Sim, struct data *Data, double *theta, 
       int nt, int nr, int *status)
{
  struct data_individual *Ind = &Data->Treatment[nt].Individual[nr];
  int n = nsteps(Sim, Ind);
  // int *t = integervector(Ind->np);
  // double **y = doublematrix(Ind->np, 1+Sim->V.variables);
  double **y = doublematrix(n, 1+Sim->V.variables);
  double logL;
  TimePoints T;

  /* ------------------------------------------------------------ setup timepoints */
  T.mode = FIT;
  T.n = Ind->np;
  T.i = 0;
  T.t = Ind->Y;
  // T.t = doublevector(T.n);
  // for (i = 0; i < T.n; i++)
  //   T.t[i] = Ind->Y[i][0];
  if ((*status = Model(nt, nr, n, theta, y, &T)) == SUCCESS) {
    // logL = 0;
    // for (i = 0; i < Ind->np; i++)
      // logL += logLikelihood(nt, nr, theta, y[t[i]], Ind->Y[i], Ind->value[i]);
      logL = logLikelihood(nt, nr, theta, y, &T);
  }
  else
    logL = -INFINITY;

  // free(t); 
  free(y[0]); free(y); 
  // free(T.t);
  return logL;
}

/* ------- return the value of the sum of a set of hyperparameters of treatment nt */
double Hypervalue(Hyperblock *HB, double *hyper, int pindex, int nt)
{
  int i, j;
  double value = 0;

  j = In_set(pindex, HB->associated_par, HB->nassociated_pars);
  if (j == -1) {
    fprintf(stderr, "parameter %d, treatment %d combination not found", pindex, nt);
    Error("");
  }
  for (i = 0; i < HB->nhypers_with_treatment[j][nt]; i++)
    value += hyper[HB->hyper_with_treatment[j][nt][i]];
  return value;
}

double sumlogPriors(Simulation *Sim, int nt, double *p, double *logP,
        struct parameterset *Q, double *hyper)
{
  int k;
  double sumlogP = 0;

  for (k = 0; k < Q->parameters; k++) 
    if (Use(FITTED, Q->P[k].type))
      sumlogP += logP[k] = logPrior(Sim, nt, p[k], &Q->P[k], hyper);

  return sumlogP;
}

double logPrior(Simulation *Sim, int nt, double p, struct parameter *P, double *hyper)
{
  double shape, location, scale, pp;
  Hyperblock *HB = Sim->HP.block;

  if (P->number == REAL)
    pp = p;
  else 
    pp = rint(p);

  if (P->hashyper) {
    if (P->distribution == NORMAL || P->distribution == TRUNCNORMAL) {
      location = Hypervalue(&HB[P->hyper_location_block_index], hyper, P->index, nt);
      scale    = Hypervalue(&HB[P->hyper_scale_block_index], hyper, P->index, nt);
    }
    else
      Error("logprior");
  }
  else if (P->distribution != USERDEFINED) {
    shape = P->shape;
    scale = P->scale;
    location = P->location;
  }

  if (P->distribution == NORMAL || P->distribution == NORMALMEAN)
    return normalPDF(pp, location, scale);
  else if (P->distribution == TRUNCNORMAL)
    return trnormalPDF(pp, location, scale);
  else if (P->distribution == UNIFORM)
    return uniformPDF(pp, location, scale);
  else if (P->distribution == EXPONENTIAL)
    return exponentialPDF(pp, location);
  else if (P->distribution == GAMMA)
    return gamma_scalePDF(pp, shape, scale);
  else if (P->distribution == BETA)
    return betaPDF(pp, shape, scale);
  else if (P->distribution == BINOMIAL)
    return binomialPMF(pp, shape, scale);
  else if (P->distribution == POISSON)
    return poissonPMF(pp, location);
  else if (P->distribution == UNIFORM0)
    return uniformPDF(pp, 0, scale);
  else if (P->distribution == GEOMETRIC)
    return geometricPDF(pp, location);
  else if (P->distribution == USERDEFINED)
    return UserPDF(pp);
  else {
    printf("%d %d\n", P->distribution, NORMAL);
    Error("distribution not implemented");
  }
  return 0;
}

void Proposal(gsl_rng *stream, double *a, struct parameter_block *PB)
{
  int i, j, npar = PB->npar, *pblock = PB->pblock;
  double z[npar], da;

  for (i = 0; i < npar; i++)
    z[i] = gsl_ran_gaussian(stream, 1.);

  for (i = 0; i < npar; i++) {
    da = 0;
    for (j = 0; j <= i; j++) 
      da += gsl_matrix_get(PB->L, i, j)*z[j];

    a[pblock[i]] += da;
  }
}

void AdjustCovariance(Simulation *Sim, struct parameter_block *PB, 
          struct parameterset *Q, int phase) 
{
  double rate;

  rate = PB->accept/PB->proposed;

  if (rate == 0 && PB->covstart < Sim->mincovsamples)
    /* -------------------------------------- scale down if not collecting samples */
    gsl_matrix_scale(PB->L, 1./Sim->scalefactor);
  else if (PB->covstart >= Sim->mincovsamples) {
    /* -------------------- rescale covariance matrix according to acceptance rate */
    if (rate < Sim->lowrate) {
      PB->covscale /= Sim->scalefactor;
      PB->inrange = 0;
      if (PB->covscale < 1./Sim->covthresh)  {
  fprintf(stderr, "%e %e\n", PB->covscale, 1./Sim->covthresh);
  Error("Scale of covariance matrix too small. See manual for causes");
      }
    }
    else if (rate > Sim->highrate) {
      PB->covscale *= Sim->scalefactor;
      PB->inrange = 0;
      if (PB->covscale > Sim->covthresh) {
  if (phase == MINTEMP) {
    PB->error = TRUE;
    PB->done = TRUE;
  }
  else {
    fprintf(stderr, "Scale of covariance matrix too large. Try reducing"
      " initscale to %f in parameter file", Sim->initscale/10.);
    Error("");
  }
      }
    }
    else
      PB->inrange++;

    /* -------------------------------- start storeing for final covariance matrix */
    if (!PB->covfinal && PB->inrange >= Sim->range) {
      PB->covfinal = TRUE;
      PB->storeno = PB->covstart;
    }
    CovMatrix(Sim, PB);
  }
}

void CovMatrix(Simulation *Sim, struct parameter_block *PB)
{
  int i, j, k, start, end;
  int  npar = PB->npar;
  double m[npar], x;
  gsl_matrix *tmp = gsl_matrix_alloc(npar, npar);
  gsl_matrix *L = PB->L;

  /* ------------------------------------------ temporarily save covariance matrix */
  gsl_matrix_memcpy(tmp, L);

  /* ------------------------------------------------- calculate covariance matrix */
  /* ------------------------------------------------------------ lower triangular */
  if (PB->covfinal) {
    start = MAX(0, PB->storeno-Sim->covsamples);
    end = PB->storeno;
  } else {
    start = 0;
    end = PB->covstart;
  }

  for (i = 0; i < npar; i++)
    m[i] = mean(PB->theta_store[i], start, end);

  for (j = 0; j < npar; j++)
    for (i = j; i < npar; i++) {
      x = 0;
      for (k = start; k < end; k++)
        x += (PB->theta_store[i][k]-m[i])*(PB->theta_store[j][k]-m[j]);
      x *= PB->covscale/(double)(end-start);
      gsl_matrix_set(L, i, j, x);
    }

  /* ------------------------------------------------------ Cholesky decomposition */

  if (gsl_linalg_cholesky_decomp(L) == GSL_EDOM) {
    printf("CHOLESKY ");
    /* -------------------------------------------- restore matrix if did not work */
    gsl_matrix_memcpy(L, tmp);
  }

  gsl_matrix_free(tmp);
}

void ClearThread(struct thread *Thread, Simulation *Sim, struct data *Data, 
     int ntp, int nrp)
{
  int c, nt, nr;
  struct sim_treatment *TT;

  if (Sim->HP.nhyper_parameters)
    free(Thread->hyper);

  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
    TT = &Thread->Treatment[nt];
    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      for (c = 0; c < Sim->Treatment[nt].Individual[nr].chains; c++) {
  free(TT->Individual[nr].Chain[c].theta);
  free(TT->Individual[nr].Chain[c].logP);
      }
    free(TT->Individual[nr].Chain);
    }
    free(TT->Individual);
  }
  free(Thread->Treatment);

  gsl_rng_free(Thread->stream);
}

void SetupThread(struct thread *Thread, Simulation *Sim, struct data *Data, 
     int phase, int ntp, int nrp, int thread)
{
  int i, c, k, nt, nr, np;
  double logL;
  struct sim_individual *Ind;
  struct chain *Chain, *Chain0, *ChainT;
  struct sim_treatment *T, *TT;

  /* ---------------------------------------------------- initialise thread hypers */
  if (Sim->HP.nhyper_parameters)
    Thread->hyper = doublevector(Sim->HP.nhyper_parameters);

  /* ---------------------------------------------------- initialise thread chains */
  Thread->Treatment = (struct sim_treatment*)
    malloc(Data->ntrt*sizeof(struct sim_treatment));

  Thread->stream = gsl_rng_alloc(gsl_rng_taus);    
  gsl_rng_set(Thread->stream, Sim->seed+thread);
  
  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
    T = &Sim->Treatment[nt];
    TT = &Thread->Treatment[nt];
    TT->Individual = (struct sim_individual*) 
      malloc(Data->Treatment[nt].nind*sizeof(struct sim_individual));
    

    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      Ind = &TT->Individual[nr];
      Ind->chains = Sim->Treatment[nt].Individual[nr].chains;
      Ind->Chain = (struct chain*) malloc(Ind->chains*sizeof(struct chain));

      np = Sim->Treatment[nt].Individual[nr].Q.parameters;
      for (c = 0; c < Ind->chains; c++) {
  Ind->Chain[c].theta = doublevector(np);
  Ind->Chain[c].logP = doublevector(np);
      }
    }
  }


  for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt]) {
    T = &Sim->Treatment[nt];
    TT = &Thread->Treatment[nt];

    for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
      Ind = &T->Individual[nr];
      Chain0 = &Ind->Chain[0];

      /* ----------------------------------------------- copy parameters from mode */
      if (phase == MINTEMP || phase == ADAPTIVE_INDIVIDUAL)
  for (c = 0; c < Ind->chains; c++) {
    Chain = &Ind->Chain[c];
    ChainT = &TT->Individual[nr].Chain[c]; 

      ChainT->logL = Chain->logL; 
      ChainT->logpost = Chain->logpost; 
      for (k = 0; k < Ind->Q.parameters; k++) {
        ChainT->theta[k] = Chain->theta[k]; 
        if (Use(FITTED, Ind->Q.P[k].type)) ChainT->logP[k] = Chain->logP[k]; 
      }
    } 
      /* -------------------------------------- copy parameters from random sample */
      else {
  if (Chain0->loop)
    i = uniform_intRNG(0, MAXSTORE, Thread->stream);
  else {
    i = uniform_intRNG(0, Chain0->laststoreno, Thread->stream);
  }
  logL = Chain0->logLstore[i];
    
  for (c = 0; c < Ind->chains; c++) {
    Chain = &Ind->Chain[c];
    ChainT = &TT->Individual[nr].Chain[c]; 
      
    for (k = 0; k < Ind->Q.parameters; k++) 
      if (Use(FITTED_AUXILLARY, Ind->Q.P[k].type))
        ChainT->theta[k] = Chain0->theta_store[k][i];
      else if (Use(CONSTANT, Ind->Q.P[k].type))
        ChainT->theta[k] = Chain->theta[k];
      else if (Use(INDIVIDUAL_CONSTANT, Ind->Q.P[k].type))
        ChainT->theta[k] = Chain->theta[k];
    
    ChainT->logL = logL;
    ChainT->logpost = Ind->temp[c]*logL + 
      sumlogPriors(Sim, nt, ChainT->theta, ChainT->logP, &Ind->Q, Sim->HP.hyper);
  }
      }
    }
  }

  /* ------------------------------------------ copy hyper parameters into threads */
  for (k = 0; k < Sim->HP.nhyper_parameters; k++) if (Sim->HP.H[k].use)
    Thread->hyper[k] = Sim->HP.hyper[k];
}

void SaturatedLikelihood(Simulation *Sim, struct data *Data, int nvar)
{
  int nt, nr;
  // double *Y;
  double **Y;
  struct data_individual *D;
  // struct chain *C;
  // TimePoints T;


  if (Sim->verbose) printf("Saturated likelihoods\n");
  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      D = &Data->Treatment[nt].Individual[nr];

      // T.n = D->np;
      // T.t = D->Y;

      if (finite(Sim->saturatedlogL))
        D->saturated = Sim->saturatedlogL;
      else {
        // Y = doublevector(nvar);
        Y = doublematrix(D->np, nvar);
        // C = &Sim->Treatment[nt].Individual[nr].Chain[0];

        D->saturated = 0;
      // for (i = 0; i < D->np; i++) {
        // SaturatedModel(nt, nr, Y, C->theta, &T);
        // D->saturated += logLikelihood(nt, nr, C->theta, Y, &T);
        // SaturatedModel(nt, nr, Y, C->theta, D->Y[i], D->value[i]);
        // D->saturated += logLikelihood(nt, nr, C->theta, Y, D->Y[i], D->value[i]);
      // }
        free(Y[0]);
        free(Y);
      }
      if (Sim->verbose) printf("%d %d: %.2f\n", nt, nr, D->saturated);
    }
}

int Time_Index(struct data_individual *Ind, int nsteps, double **y, int *t_index)
{
  int j, t = 0;

  for (j = 0; j < Ind->np; j++) {
    while(t < nsteps && y[t][0] < Ind->Y[j][0]) t++;
    if (t == nsteps) {
      printf("index model\n");
      for (j = 0; j < nsteps; j++)
  printf("%d %f\n", j, y[j][0]);
      printf("index data\n");
      for (j = 0; j < Ind->np; j++) 
  printf("%d %f\n", j, Ind->Y[j][0]);
      Error("data and model times do not match");
    }
    if (y[t][0] == Ind->Y[j][0])
      t_index[j] = t;
    else if (t == 0) {
      printf("index model\n");
      for (j = 0; j < nsteps; j++)
  printf("%d %f\n", j, y[j][0]);
      printf("index data\n");
      for (j = 0; j < Ind->np; j++) 
  printf("%d %f\n", j, Ind->Y[j][0]);
      Error("data and model times do not match");
      return FALSE;
    }
    else if (Ind->Y[j][0]-y[t-1][0] < y[t][0]-Ind->Y[j][0])
      t_index[j] = t-1;
    else
      t_index[j] = t;
  }
  return TRUE;
}

int nsteps(Simulation *Sim, struct data_individual *Ind)
{
  double start, end;

  if (Sim->notime) return Ind->np;
  if (isinf(Sim->timestart))
    start = Ind->Y[0][0];
  else
    start = Sim->timestart;
  if (isinf(Sim->timeend))
    end = Ind->Y[Ind->np-1][0];
  else
    end = Sim->timeend;
  return lrint((end-start)/timestep())+1;
}

void Best(Simulation *Sim, struct data *Data, struct parameterset *Q, 
    struct chain *ChainT, struct chain *Chain,
    int ntp, int nrp, int c, int phase)
{
  int k;

  Chain->logpost = ChainT->logpost;
  Chain->logL = ChainT->logL;
  for (k = 0; k < Q->parameters; k++) if (Use(FITTED, Q->P[k].type))
    Chain->theta[k] = ChainT->theta[k];
}

int Done(Simulation *Sim, struct data *Data, int ntp, int nrp, int phase)
{
  int nt, nr, c, b;
  struct sim_individual *Ind;

  if (phase == MINTEMP || phase == ADAPTIVE_INDIVIDUAL) {
    for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt])
      for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr]) {
  Ind = &Sim->Treatment[nt].Individual[nr];
  for (c = 0; c < Ind->chains; c++)
    for (b = 0; b < Ind->Q.Blocks.nblocks; b++)
      if (!Ind->Chain[c].PB[b].categorical)
        if (!Ind->Chain[c].PB[b].done) return FALSE;
      }
  }
  else if (phase == NONADAPTIVE)
    for (nt = NTS; nt <= NTE; nt++) if (Sim->trtused[nt])
      for (nr = NRS; nr <= NRE; nr++) if (Sim->indused[nt][nr])
  for (c = 0; c < Sim->Treatment[nt].Individual[nr].chains; c++)
    if (!Sim->Treatment[nt].Individual[nr].Chain[c].done) return FALSE;
  
  return TRUE;
}

void SetNAN(double *x, int n)
{
  int i;
  for (i = 0; i < n; i++)
    x[i] = NAN;
}

void Default(Simulation *Sim)
{
  nt1 = nr1            = 0;
  nt2 = nr2            = -1;
  Sim->rdp             = FALSE;
  Sim->sim             = FALSE;
  Sim->thin            = 100;
  Sim->mcmc            = TRUE;
  Sim->popd            = TRUE;
  Sim->marg            = TRUE;
  Sim->seed            = 0;
  Sim->xrate           = 0;
  Sim->range           = 2;
  Sim->error           = 5;
  Sim->notime          = FALSE;
  Sim->chains          = 1;
  Sim->copyind         = -1;
  Sim->threads         = omp_get_num_procs();
  Sim->verbose         = TRUE;
  Sim->timeend         = INFINITY;
  Sim->samples         = 1000;
  Sim->derived         = 0;
  Sim->discard         = 1000;
  Sim->lowrate         = 0.2;
  Sim->highrate        = 0.4;
  Sim->minprior        = -3.;
  Sim->initscale       = 500.;
  Sim->covthresh       = 1e30;
  Sim->gsl_error       = FALSE;
  Sim->timestart       = -INFINITY;
  Sim->factors[0]      = '\0';
  Sim->covsamples      = 500;
  Sim->population      = 0;
  Sim->growthcurves    = FALSE;
  Sim->scalefactor     = 2.;
  Sim->parsperblock    = 50;
  Sim->phase1adjust    = 100;
  Sim->phase2adjust    = 1000;
  Sim->mincovsamples   = 50;
  Sim->ignorewarning   = FALSE;
  Sim->usehypermeans   = TRUE;
  Sim->minexchangerate = 0.0;
  Sim->saturatedlogL   = INFINITY;
}

double GelmanRubinStats(Simulation *Sim, struct data *Data, FILE *fp)
{
  int nt, nr, c, i, k, n, count, *cc = integervector(Sim->threads);
  double B, W, R, psi, *psi_c = doublevector(Sim->threads);
  double *s2 = doublevector(Sim->threads), maxR = 0;
  struct sim_individual *Ind;
  struct chain *Chain;

  if (Sim->threads == 1) return 0;

  for (nt = nt1; nt <= NT2; nt++) if (Sim->trtused[nt])
    for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
      Ind = &Sim->Treatment[nt].Individual[nr];
      Chain = &Ind->Chain[0];

      for (c = 0; c < Sim->threads; c++)
  cc[c] = 0;
      for (i = 0; i < Chain->storeno; i++)
  cc[Ind->chainno[i]]++;
      n = Chain->storeno;
      for (c = 0; c < Sim->threads; c++)
  if (cc[c] < n) n = cc[c];
      
      for (k = 0; k < Ind->Q.parameters; k++) 
  if (Use(FITTED, Ind->Q.P[k].type)) {
    psi = B = W = 0;
    for (c = 0; c < Sim->threads; c++)
      psi_c[c] = s2[c] = 0;
    
    for (c = 0; c < Sim->threads; c++) {
      count = 0;
      for (i = 0; i < Chain->storeno && count < n; i++)
        if (Ind->chainno[i] == c) {
    psi_c[c] += Chain->theta_store[k][i];
    count++;
        }
    }
    for (c = 0; c < Sim->threads; c++) {
      psi_c[c] /= (double)n;
      psi += psi_c[c];
    }
    psi /= (double)Sim->threads;
    for (c = 0; c < Sim->threads; c++)
      B += (psi_c[c]-psi)*(psi_c[c]-psi);
    
    for (c = 0; c < Sim->threads; c++) {
      count = 0;
      for (i = 0; i < Chain->storeno && count < n; i++)
        if (Ind->chainno[i] == c) {
    s2[c] += (Chain->theta_store[k][i]-psi_c[c])*
      (Chain->theta_store[k][i]-psi_c[c]);
    count++;
        }
    }
    for (c = 0; c < Sim->threads; c++)
      W += s2[c];
    
    R = sqrt((n-1.)*(1./n+Sim->threads*B/((Sim->threads-1.)*W)));
    if (R > maxR) maxR = R;
    if (fp) fprintf(fp, "%d %d %d %f\n", nt, nr, k, R);
  }
    }
  if (fp && maxR > 1.1) fprintf(stderr, "WARNING! Poor convergence of Markov chains. "
        "Use 'bag j' to examine Gelman Rubin statistics.\n");

  free(s2); free(psi_c); free(cc);
  return maxR;
}

void ReadParameters(Simulation *Sim, struct parameterset *Q, 
        int nt, int nr, double *p)
{
  int n, i, j, k, countpars = 0, countzero = 0;
  double v;
  FILE *fp = fopen(Sim->parsfile, "r");
  static int firsttime = TRUE;
  struct parameter *P = Q->P;
  struct hyper *H;
  double *hyper = Sim->HP.hyper;
  Hyperblock *HB = Sim->HP.block;

  /* ------------------------------------------------------------ hyper parameters */
  /* ---------------------------------------------- take location as initial value */
  /* --------------------------- must do this before setting individual parameters */
  if (Sim->HP.nhyper_parameters) {
    for (i = 0; i < Sim->HP.nhyper_parameters; i++) if (Sim->HP.H[i].use) {
      H = &Sim->HP.H[i];
      if (H->hdist == GAMMA || H->hdist == INVGAMMA || H->hdist == PARETO)
  hyper[i] = H->scale;
      else if (H->hdist == BETA) 
  hyper[i] = H->shape/(H->shape+H->scale);
      else
  hyper[i] = H->location;
    }
  }

  /* ------------------------------------------------------- individual parameters */
  /* ----------------------------------- first, use location of prior distribution */
  for (k = 0; k < Q->parameters; k++) {
    if (Use(CONSTANT, P[k].type))
      p[k] = P[k].location;
    else if (P[k].hashyper)
      if (P[k].distribution == GAMMA || P[k].distribution == UNIFORM0)
  p[k] = Hypervalue(&HB[P[k].hyper_scale_block_index], hyper, P[k].index, nt);
      else if (P[k].distribution == BINOMIAL)
  p[k] = P[k].shape*Hypervalue(&HB[P[k].hyper_scale_block_index], hyper, 
             P[k].index, nt);
      else
  p[k] = Hypervalue(&HB[P[k].hyper_location_block_index], hyper, P[k].index,
        nt);
    else if (Use(FITTED, P[k].type)) {
      if (P[k].distribution == UNIFORM)
  p[k] = 0.5*(P[k].location+P[k].scale);
      else if (P[k].distribution == UNIFORM0)
  p[k] = 0.5*P[k].scale;
      else if (P[k].distribution == GAMMA || P[k].distribution == INVGAMMA)
  p[k] = P[k].shape*P[k].scale;
      else if (P[k].distribution == BETA)
  p[k] = P[k].shape/(P[k].shape+P[k].scale);
      else if (P[k].distribution == BINOMIAL)
  p[k] = P[k].shape*P[k].scale;
      else if (P[k].distribution == CATEGORICAL) {
  v = 0;
  for (i = 0; i < P[k].ncategories; i++) 
    if (P[k].p[i] > v) {
      v = P[k].p[i];
      j = i;
    }
  p[k] = j;
      }
      else if (P[k].distribution == USERDEFINED)
  p[k] = 0;
      else
  p[k] = P[k].location;
    }
}
  /* ------------------------------------ second, read from best file if it exists */

  if (fp) {
    rewind(fp);
    do {
      n = fscanf(fp, "%d%d%d%lf", &i, &j, &k, &v);
      if (n == EOF) break;
      if (n != 4) {
  fprintf(stderr, "%s is in wrong format: treatment individual parameter "
    "value\n", Sim->parsfile);
  Error("");
      }
      if (k >= Q->parameters) {
  if (firsttime)
    fprintf(stderr, "parameter %d not read in for treatment %d, individual %d "
      "as it is not used\n", k, nt, nr);
      }
      else { 
  if (!Use(CONSTANT, P[k].type) && P[k].type != NOTUSED && i == nt &&
      ((Sim->copyind == -1 && j == nr) || Sim->copyind == j))
    p[k] = v;
  if (Use(INDIVIDUAL_CONSTANT, P[k].type)) 
    P[k].location = v;
      }
    } while (TRUE);
    rewind(fp);
  }

  /* -------------------------------------- if parameter has a hyper use hypermean */
  /* if (Sim->usehypermeans) */
  /*   for (k = 0; k < Q->parameters; k++) { */
  /*     if (P[k].hashyper) { */
  /*   if ((i=In_set(k, Sim->H.par_indices, Sim->H.npar)) != -1) */
  /*     p[k] = Sim->H.hypermean.a[i]; */
  /*   else if ((i=In_set(k, T->H.par_indices, T->H.npar)) != -1) */
  /*     p[k] = T->H.hypermean.a[i]; */
  /*   else */
  /*     Error("hhh"); */
  /*   //  printf("%d %d %d %e\n", nt, nr, k, p[k]); */
  /*     } */
  /*   } */

  /* ------------------------------------------------- third, read from model file */
  Parameters(nt, nr, p, Q);
  for (k = 0; k < Q->parameters; k++) 
    if (Use(INDIVIDUAL_CONSTANT, P[k].type)) 
      P[k].location = p[k];

  /* ----------------------------------------- check that parameters have been set */
  for (k = 0; k < Q->parameters; k++) 
    if (P[k].distribution != CATEGORICAL && Use(FITTED, P[k].type)) {
      countpars++;
      if (p[k] == 0) countzero++;
    }
  if (!Sim->ignorewarning && countpars > 4 && countzero > countpars/4) {
    fprintf(stderr, "many parameters are 0 for treatment %d individual %d\n", 
      nt, nr);
    for (k = 0; k < Q->parameters; k++) if (Use(FITTED, P[k].type))
      fprintf(stderr, "%d %e\n", k, p[k]);
    Error("");
  }

  /* ----------------------------------- check that parameters have correct values */
  for (k = 0; k < Q->parameters; k++) if (Use(FITTED, P[k].type)) {
    if (isinf(p[k]) || isnan(p[k])) {
      fprintf(stderr, "parameter %d=%e for treatment %d, individual %d\n", 
        k, p[k], nt, nr);
      Error("");
    }
  }

  if (fp) fclose(fp);
  firsttime = FALSE;
}

/* -------------------------------------------------------- clear solution vectors */
void Clear(double **y, int n, int m)
{
  int i,  j;
  for (i = 0; i < n; i++)
    for (j = 1; j < m; j++)
      y[i][j] = 0;
  for (i = 0; i < n; i++)
    y[i][0] = i;
}

int Use(int t1, int t2) 
{
  if ((t1 == ALL                      && t2 != NOTUSED)                     ||
      (t1 == CONSTANT                 && t2 == CONSTANT)                    ||
      (t1 == INDIVIDUAL_CONSTANT      && t2 == INDIVIDUAL_CONSTANT)         ||
      (t1 == DERIVED                  && t2 == DERIVED)                     ||
      (t1 == AUXILLARY                && t2 == AUXILLARY)                   ||
      (t1 == FITTED                   && t2 == FITTED)                      ||
      (t1 == FITTED_CONSTANT          && 
       (t2 == FITTED || t2 == CONSTANT || t2 == INDIVIDUAL_CONSTANT))       ||
      (t1 == FITTED_DERIVED           && (t2 == FITTED || t2 == DERIVED))   ||
      (t1 == FITTED_AUXILLARY         && (t2 == FITTED || t2 == AUXILLARY)) ||
      (t1 == FITTED_DERIVED_AUXILLARY &&
       (t2 == FITTED || t2 == DERIVED || t2 == AUXILLARY)))
    return TRUE;
  else
    return FALSE;
}

void Gibbs(Simulation *Sim, struct data *Data, struct thread *T, gsl_rng *stream)
{
  int i, b, k, j, n, nt, nr, hyperdist, paramdist, hp, hp1, j1, p;
  int nhp = Sim->HP.nhyper_parameters;
  double sum, s2, mu, r;
  double newvalue[nhp];
  struct hyper *H = Sim->HP.H;
  Hyperparameters *HP = &Sim->HP;
  Hyperblock HB;
  double tau[nhp], mean[nhp], shape[nhp], scale[nhp];

  for (b = 0; b < HP->nblocks; b++) if (HP->block[b].use) {
    HB = HP->block[b];
    n = HB.nhypers_used;
    paramdist = HB.pdistribution;
    hyperdist = HB.hdistribution;

    if (paramdist == NORMAL || paramdist == TRUNCNORMAL || paramdist == NORMALMEAN) {
      if (hyperdist == NORMAL) {
  /* ------------------------------------------------------ mean of a normal */
  /* ------------------------ set up tau and mean for hypers used in block i */
  // tau = Sigma_0^-1 + n*Sigma^-1
  // mean = tau^-1 * (Sigma_0^{-1}*mu_0 + n*Sigma^{-1}*xbar)
  for (j = 0; j < n; j++) {
    hp = HB.hyper_used[j];
    mean[hp] = H[hp].location/H[hp].scale;
    tau[hp] = 1./H[hp].scale;
  }

  /* ------------------- loop through all parameters associated with block i */
  for (j = 0; j < HB.nassociated_pars; j++) {
    p = HB.associated_par[j];
    
    /* -- loop through each treatment associated with parameter p in block i */
    for (k = 0; k < HB.ntreatments_with_parameter[j]; k++) {
      nt = HB.treatment_with_parameter[j][k];
      if (!Sim->trtused[nt]) continue;

      /* ---------------- variance hyperparameter for treatment nt of block i*/
      s2 = Hypervalue(&HP->block[HB.dependent_block[j]], T->hyper, p, nt);

      /* ------------ sum the parameter p over all individuals in treatment nt */
      sum = 0;
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr])
        sum += T->Treatment[nt].Individual[nr].Chain[0].theta[p];

      /* --- in each treatment, loop through all hyperparameter combinations */
      for (j1 = 0; j1 < HB.nhypers_with_treatment[j][nt]; j1++) {
        hp1 = HB.hyper_with_treatment[j][nt][j1];
        mean[hp1] += sum/s2;
        tau[hp1] += HP->nind[nt]/s2;
      }
    }
  }

  for (j = 0; j < n; j++) {
    hp = HB.hyper_used[j];
    newvalue[hp] = normalRNG(mean[hp]/tau[hp], 1./tau[hp], stream);
  }
      }
      else if (hyperdist == GAMMA || hyperdist == INVGAMMA) {
  /* ----------- initialise shape and scale for each variance hyperparameter */
  for (j = 0; j < n; j++) {
    hp = HB.hyper_used[j];
    shape[hp] = 2.*H[hp].shape;
    scale[hp] = 2.*H[hp].scale;
  }
  
  /* ------------------- loop through all parameters associated with block i */
  for (j = 0; j < HB.nassociated_pars; j++) {
    p = HB.associated_par[j];
    
    /* -- loop through each treatment associated with parameter p in block i */
    for (k = 0; k < HB.ntreatments_with_parameter[j]; k++) {
      nt = HB.treatment_with_parameter[j][k];
      if (!Sim->trtused[nt]) continue;
      
      mu = Hypervalue(&HP->block[HB.dependent_block[j]], T->hyper, p, nt);
      sum = 0;
      for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) {
        r = T->Treatment[nt].Individual[nr].Chain[0].theta[p]-mu;
        sum += r*r;
      }

      /* ---------------- in each treatment, loop through all hyperparameter */
      for (i = 0; i < HB.nhypers_with_treatment[j][nt]; i++) {
        hp = HB.hyper_with_treatment[j][nt][i];
        shape[hp] += HP->nind[nt];
        scale[hp] += sum;
       }
    }
  }

  /* ------------------------------------------- randomly draw new variances */
  for (j = 0; j < n; j++) {
    hp = HB.hyper_used[j];
    if (hyperdist == GAMMA)
      newvalue[hp] = gammaRNG(shape[hp]/2., scale[hp]/2., stream);
    else
      newvalue[hp] = invgammaRNG(shape[hp]/2., scale[hp]/2., stream);
  }
      }
    }
    /* else if (paramdist == EXPONENTIAL) { */
    /*   /\* -------------------------------------------------- rate of an exponential *\/ */
    /*   shape = doublevector(nhp); */
    /*   scale = doublevector(nhp); */

    /*   /\* ---------------------- initialise shape and scale for each hyperparameter *\/ */
    /*   for (j = 0; j < n; j++) { */
    /*   hp = HB.hyper_used[j]; */
    /*   shape[hp] = H[hp].shape; */
    /*   scale[hp] = H[hp].scale; */
    /*   } */
  
    /*   /\* --------------------- loop through all parameters associated with block i *\/ */
    /*   for (j = 0; j < HB.nassociated_pars; j++) { */
    /*   p = HB.associated_par[j]; */
    
    /*   /\* ---- loop through each treatment associated with parameter p in block i *\/ */
    /*   for (k = 0; k < HB.ntreatments_with_parameter[j]; k++) { */
    /*     nt = HB.treatment_with_parameter[j][k]; */
    /*     if (!Sim->trtused[nt]) continue; */
     
    /*     /\* ----------------- in each treatment, loop through all hyperparameters *\/ */
    /*     for (i = 0; i < HB.nhypers_with_treatment[j][nt]; i++) { */
    /*       hp = HB.hyper_with_treatment[j][nt][i]; */
    /*       shape[hp] += HP->nind[nt]; */
    /*       for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) */
    /*         scale[hp] += T->Treatment[nt].Individual[nr].Chain[0].theta[p]; */
    /*     } */
    /*   } */
    /*   } */
    /*   /\* ------------------------------------------------- randomly draw new rates *\/ */
    /*   for (j = 0; j < n; j++) { */
    /*   hp = HB.hyper_used[j]; */
    /*   newvalue[hp] = gammaRNG(shape[hp], scale[hp], stream); */
    /*   } */
    /*   free(scale); free(shape); */
    /* } */
    /* else if (paramdist == POISSON) { */
    /*   /\* ------------------------------------------------------ scale of a poisson *\/ */
    /*   shape = doublevector(nhp); */
    /*   N = integervector(nhp); */

    /*   /\* -------------------------- initialise shape and N for each hyperparameter *\/ */
    /*   for (j = 0; j < n; j++) { */
    /*   hp = HB.hyper_used[j]; */
    /*   shape[hp] = H[hp].shape; */
    /*   N[hp] = 0; */
    /*   } */
  
    /*   /\* --------------------- loop through all parameters associated with block i *\/ */
    /*   for (j = 0; j < HB.nassociated_pars; j++) { */
    /*   p = HB.associated_par[j]; */
    
    /*   /\* ---- loop through each treatment associated with parameter p in block i *\/ */
    /*   for (k = 0; k < HB.ntreatments_with_parameter[j]; k++) { */
    /*     nt = HB.treatment_with_parameter[j][k]; */
    /*     if (!Sim->trtused[nt]) continue; */
     
    /*     /\* ----------------- in each treatment, loop through all hyperparameters *\/ */
    /*     for (i = 0; i < HB.nhypers_with_treatment[j][nt]; i++) { */
    /*       hp = HB.hyper_with_treatment[j][nt][i]; */
    /*       N[hp] += HP->nind[nt]; */
    /*       for (nr = nr1; nr <= NR2; nr++) if (Sim->indused[nt][nr]) */
    /*         shape[hp] += T->Treatment[nt].Individual[nr].Chain[0].theta[p]; */
    /*     } */
    /*   } */
    /*   } */
    /*   /\* ------------------------------------------------ randomly draw new scales *\/ */
    /*   for (j = 0; j < n; j++) { */
    /*   hp = HB.hyper_used[j]; */
    /*   newvalue[hp] = gammaRNG(shape[hp], H[hp].scale/(N[hp]*H[hp].scale+1.), */
    /*         stream); */
    /*   } */
    /*   free(shape); free(N); */
    /* } */
  }
  for (hp = 0; hp < nhp; hp++)
    if (H[hp].use) T->hyper[hp] = newvalue[hp];
}

void ConstructPiecewiseLine(struct line *L)
{
  int c;

  for (c = 0; c < L->n-1; c++) {
    L->m[c] = (L->y[c]-L->y[c+1])/(L->x[c]-L->x[c+1]);
    L->c[c] = L->y[c]-L->m[c]*L->x[c];
  }
}

void ReturnPiecewiseLine(double x, double *f, double *df, struct line *L)
{
  int c = 0;

  if (L->x[0] > L->x[1]) 
    while (x < L->x[c+1]) c++;
  else
    while (x > L->x[c+1]) c++;
  *f = L->m[c]*x+L->c[c];
  *df = L->m[c];
}

void Used(int *used, char *line) 
{
  const char *ptr = line;
  char field[1000];
  int n, j, seq = 0, last = -1, i;

  while (sscanf(ptr, "%999[^-,]%n", field, &n) == 1) {
    if (strcmp("all", field) == 0) {
      for (i = 0; i < MAXTRTSREPS; i++)
  used[i] = TRUE;
      break;
    }
    j = atoi(field);
    if (j >= MAXTRTSREPS) {
      Error("Maximum number of treatments/individuals is MAXTRTSREPS");
    }
    ptr += n;
    if (seq) {
      for (i = last; i <= j; i++)
  used[i] = TRUE;
      seq = 0;
    }
    if (*ptr == '-') {
      last = j;
      seq = 1;
    }
    else if (*ptr == ',' || *ptr == '\0')
      used[j] = TRUE;
    if (*ptr == '\0') break;
    ++ptr;
  }
}

void Ignore(int **used, char *line)
{
  const char *ptr = line;
  char field[1000];
  int n, j, t;

  while (sscanf(ptr, "%31[^,;]%n", field, &n) == 1) {
    j = atoi(field);
    ptr += n;
    if (*ptr == ';' || *ptr == '\0')
      used[t][j] = FALSE;
    else if (*ptr == ',')
      t = j;
    if (*ptr == '\0') break;
    ++ptr;
  }
}

int In_set(int x, int *a, int n)
{
  int i;

  if (a == NULL) return -1;

  for (i = 0; i < n; i++)
    if (x == a[i]) return i;

  return -1;
}

void SetupHyperBlock(Hyperparameters *HP, struct parameterset *Q, int k, int nt, 
         int block) 
{
  Hyperblock *HB;

  /* ----------------- check hyperparameter block is declared in HyperParameters() */
  if (block >= HP->nblocks || !HP->block[block].use) {
    fprintf(stderr, "hyper parameter block %d for treatment %d, parameter %d"
      "not declared in HyperParameters()\n", block, nt, k);
    Error("");
  }

  HB = &HP->block[block];
  if (HB->use == 1) HB->dependent_block = integervector(HB->nassociated_pars);

  /* ---------------------------- increase HB->use by 1 to indicate it is required */
  HB->use++;
  
  if (HB->pdistribution == NOTAPPLICABLE) HB->pdistribution = Q->P[k].distribution;
}
 
