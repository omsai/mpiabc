/**
   @file sampler.c

   @brief Sample parameters.

   Copyright (c) 2024, Pariksheet Nanda

   SPDX-License-Identifier: Apache-2.0
 */

#include <math.h>

#include <sampler.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>


/**
   Silverman's rule of thumb to calculate Gaussian bandwidth parameter.

   The formula used is:
   \f{eqnarray*}{
   h &=& 0.9 \min \left( \sigma,
                        {\displaystyle\frac{\mathrm{IQR}}{1.34}} \right)
            n ^{-\frac{1}{5}} \\
   \mathrm{where,}& \\
   \sigma &\equiv& \textrm{standard deviation of } x \\
   \mathrm{IQR} &\equiv& \textrm{inter-quartile range of } x \\
   n &\equiv& \textrm{cadinality of } x \\
   \f}

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
   1-D Gaussian kernel density estimate.

   The formula used is:
   \f{eqnarray*}{
   \hat{f}_h(x) &=& \frac{1}{n}
                    \sum_i \mathcal{N} ( x_i - x, h ) \\
   \mathrm{where,}& \\
   n &\equiv& \textrm{cardinality of } x \\
   \mathcal{N}(\mu, \sigma) &\equiv&
     \textrm{Normal / Gaussian distribution function} \	\
   \f}
   
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
   @param model Target model posterior density.
   @param kernel Smoothing kernel function.
   @param n Number of parameters to calibrate (maybe).
   @param sampling Sampling density.
   @param proposal Proposal density.
   @param alpha Control the effective sample size.
   @return GSL_SUCCESS

   @see kde
 */
int
abc_smc(gsl_matrix* params,
	gsl_function* model,
	gsl_function *kernel,
	int n,
	gsl_function* sampling,
	gsl_function* proposal,
	int alpha)
{
  (void)params;
  (void)model;
  (void)kernel;
  (void)n;
  (void)sampling;
  (void)proposal;
  (void)alpha;
  //gsl_sort(data->data, 1, data->size);
  return GSL_SUCCESS;
}
