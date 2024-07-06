/**
   @file sampler.c

   @brief Sample parameters.
 */

#include <math.h>

#include <sampler.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>


/**
   Silverman's rule of thumb to calculate Gaussian bandwidth parameter.

   @param x Sorted input vector.
   @return Bandwidth of Gaussian distribution.

   @see kde
*/
double
silverman(gsl_vector* x)
{
  double sd = gsl_stats_sd(x->data, 1, x->size);
  double iqr =
    gsl_stats_quantile_from_sorted_data(x->data, 1, x->size, 0.75) -
    gsl_stats_quantile_from_sorted_data(x->data, 1, x->size, 0.25);
  return 0.9 * fmin(sd, iqr / 1.34) * pow(x->size, -0.2);
}


/**
   1-D gaussian kernel density estimate.
   
   @param x Input point to evaluate the kernel density.
   @param params_ data Input vector to construct the kernel density and
                  Gaussian bandwidth parameter.
   @return Gaussian kernel density.

   @see silverman
 */
double
kde(double x, void* params_)
{
  kde_params_t* params = (kde_params_t*)params_;
  gsl_vector* data = params->data;
  double h = params->bandwidth;
  double y = 0;
  for (size_t i = 0; i < data->size; ++i) {
    y += gsl_ran_gaussian_pdf(gsl_vector_get(data, i) - x, h);
  }
  return y / data->size;
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

   @see kde
 */
int
abc_smc(gsl_block* params,
	gsl_function* target,
	gsl_function *kernel,
	int n,
	gsl_function* sampling,
	gsl_function* proposal,
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
