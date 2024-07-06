/**
   @file sampler.c

   @brief Sample parameters.
 */

#include <math.h>

#include <sampler.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>


/**
   Silverman's rule of thumb to calculate Gaussian bandwidth parameter.

   @param x Sorted input vector.
*/
double
silverman(gsl_block* x)
{
  double sd = gsl_stats_sd(x->data, 1, x->size);
  double iqr =
    gsl_stats_quantile_from_sorted_data(x->data, 1, x->size, 0.75) -
    gsl_stats_quantile_from_sorted_data(x->data, 1, x->size, 0.25);
  return 0.9 * fmin(sd, iqr / 1.34) * pow(x->size, -0.2);
}


/**
   1-D gaussian kernel density estimate.
   
   @param y Output kernel density for x.
   @param x Input vector to evaluate the kernel density for each point in data.
   @param data Input vector.
   @param h Gaussian bandwidth parameter.
   @return GSL_SUCCESS

   @see silverman
 */
int
kde(gsl_block** y,
    gsl_block* x,
    gsl_block* data,
    double h)
{
  *y = gsl_block_calloc(x->size);
  /* Ideally, one would use a sparse-FFT here, but instead the following is a
     manual convolution of our sparse data with the assumption that there are
     relatively few data points. */
  for (size_t i = 0; i < x->size; ++i) {
    for (size_t j = 0; j < data->size; ++j) {
      (*y)->data[i] += gsl_ran_gaussian_pdf(data->data[j] - x->data[i], h);
    }
  }
  for (size_t i = 0; i < x->size; ++i) {
    (*y)->data[i] /= data->size;
  }
  return GSL_SUCCESS;
}

/**
   ABC-SMC algorithm 4.8 in \cite sisson2019.

   @param params Output set of weighted parameter vectors.
   @param target Target model posterior density.
   @param kernel Smoothing kernel function.
   @param n Number of parameters to calibrate (maybe).
   @param sampling Sampling density.
   @param proposal Proposal density.
   @param alpha Control the effective sample size.
   @param sum_stat Vector of summary statistics.
   @return GSL_SUCCESS
 */
int
abc_smc(gsl_block* params,
	gsl_block* (*target)(gsl_block*),
	double (*kernel)(double),
	int n,
	double (*sampling)(double),
	double (*proposal)(double, double),
	int alpha,
	gsl_block* sum_stat)
{
  (void)params;
  (void)target;
  (void)kernel;
  (void)n;
  (void)sampling;
  (void)proposal;
  (void)alpha;
  (void)sum_stat;
  //gsl_sort(data->data, 1, data->size);
  return GSL_SUCCESS;
}
