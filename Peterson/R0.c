#include <math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

double R0(double flock_size, double housing_period, double prob_dam_to_lamb, int maxage, double *params)
{
    int n_age_classes = (maxage+1)*12;
    int c, a, i, y, t;
    int infectious_period, years_infectious, period_infectious_to_lambs;
    
    int latent_period = lrint(params[5]);
    double beta_field = params[0];
    double beta_housed = params[1];
    
    double p_mort, p_cull, sum, lambs_per_dam;
    double removal_rate[n_age_classes], infection_rate;
    double cohort_size_age_0;
    double cohort_size[maxage+1][n_age_classes];
    double infecteds[maxage+1][n_age_classes];
    double R0[n_age_classes][n_age_classes];

    // age-specific mortality and cull probabilities from Andrew's data
    double prob_mortality[] = {
        0.01,
        0.012,
        0.027,
        0.041,
        0.019,
        0.021,
        0.034,
    };

    double prob_culled[] = {
        0,
        0.046,
        0.065,
        0.093,
        0.134,
        0.196,
        0.295,
    };

    // create vector of removal rates of a cohort from ages 0 to n_age_classes
    for (a = 0; a < n_age_classes; a++) {
        p_mort = prob_mortality[a/12];
        p_cull = prob_culled[a/12];
        removal_rate[a] = -(log(1.-p_mort)+log(1.-p_cull))/12.;
    }

    // these two ways are equivalent
    // create vector of infection rates within a year starting in April
    infection_rate = ((12-housing_period)*beta_field+housing_period*beta_housed)/12.;
    // for (a = 0; a < 12; a++)
        // if (a < 12-housing_period)
        //     infection_rate[a] = beta_field;
        // else
        //     infection_rate[a] = beta_housed;

    // find intial cohort size in a non-infected flock
    // take initial cohort size as 1
    // find size of cohort as it ages
    cohort_size[0][0] = 1;
    for (a = 1; a < n_age_classes; a++)
        cohort_size[0][a] = cohort_size[0][a-1]*exp(-removal_rate[a]);
    // sum size of all existent cohorts in April
    sum = 0;
    for (y = 0; y <= maxage; y++)
        sum += cohort_size[0][12*y];
    cohort_size_age_0 = flock_size/sum;

    // for cohorts of y years old at first exposure
    // assuming that the infectious ewe is starts to be infectious in April - but
    // this is not necessary, just convenient
    for (y = 0; y <= maxage; y++) {
        // cohort size at age 0
        cohort_size[y][0] = cohort_size_age_0;

        // just removal for first y years before exposure
        for (a = 1; a < 12*y; a++)
            cohort_size[y][a] = cohort_size[y][a-1]*exp(-removal_rate[a]);
        // then removal and infection for rest of lifespan
        for (; a < n_age_classes; a++)
            cohort_size[y][a] = cohort_size[y][a-1]*exp(-removal_rate[a]-infection_rate);

        // calculate number of infecteds
        infecteds[y][0] = 0;
        for (a = 1; a < n_age_classes; a++)
            infecteds[y][a] = infection_rate/(removal_rate[a]+infection_rate)*(cohort_size[y][a-1]-cohort_size[y][a]);
    }

    // set matrix elements to zero
    for (a = 0; a < n_age_classes; a++)
        for (i = 0; i < n_age_classes; i++)
            R0[a][i] = 0;

    // for a ewe infected at age a, sum the number of ewes it infects of age i
    for (a = 0; a < n_age_classes; a++) {
        infectious_period = max(0, n_age_classes-a-latent_period);
        for (t = 0; t < infectious_period; t++) {
            // for all cohorts from past to future 
            for (y = maxage; y >= -maxage; y--) {
                // get cohort
                // +ve y are cohorts in the past
                // -ve y are future cohorts (identical to cohort 0)
                c = max(0, y);
                // age of cohort y at time t
                i = 12*y+t;
                if (0 <= i && i < n_age_classes)
                    R0[a][i] += infecteds[c][i];
            }
        }
    }

    // contribution of vertical transmission to an infectious ewe's lambs
    lambs_per_dam = cohort_size_age_0/(flock_size-cohort_size_age_0);

    for (a = 0; a < n_age_classes; a++) {
        // a ewe infected at age a is infectious for max_age-a-1-latent_period months
        infectious_period = max(0, n_age_classes-a-latent_period);
        if (a < 24)
            // only 2 year olds and older give birth
            // so ewes infectious before 24 months cannot infect their lambs 
            // until they are 2 years old
            period_infectious_to_lambs = max(0, n_age_classes-24-latent_period);
        else
            period_infectious_to_lambs = infectious_period;
        years_infectious = period_infectious_to_lambs/12;
        R0[a][0] += lambs_per_dam*prob_dam_to_lamb*years_infectious;
    }

    gsl_matrix *A = gsl_matrix_alloc(n_age_classes, n_age_classes);
    gsl_vector_complex *v = gsl_vector_complex_alloc(n_age_classes);
    gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(n_age_classes);
    for (a = 0; a < n_age_classes; a++)
        for (i = 0; i < n_age_classes; i++)
            gsl_matrix_set(A, a, i, R0[a][i]);
    gsl_eigen_nonsymm(A, v, w);
    double max_eigenvalue = 0, e;
    for (a = 0; a < n_age_classes; a++) {
        e = GSL_REAL(gsl_vector_complex_get(v, a));
        if (e > max_eigenvalue)
            max_eigenvalue = e;
    }
    gsl_eigen_nonsymm_free(w);
    gsl_vector_complex_free(v);
    gsl_matrix_free(A);

    return max_eigenvalue;
}
