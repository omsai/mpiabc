/**
   @file infer.c

   @brief Entrypoint for model inference.

   Copyright (c) 2024, Pariksheet Nanda

   SPDX-License-Identifier: Apache-2.0
 */

#include <abc.h>

#include <gsl/gsl_randist.h>


/**
   Sample from the model using ABC-SMC.
 */
int
main() {
  double a_hyperparams[] = {0, 1};
  double b_hyperparams[] = {0, 1};
  double c[] = {1.50};
  double d[] = {0.75};
  gsl_ran_function gaussian_tail;
  gaussian_tail.func2p = gsl_ran_gaussian_tail;
  param_t params[] = {
    {"a", &gaussian_tail, a_hyperparams, 2},
    {"b", &gaussian_tail, b_hyperparams, 2},
    {"c", NULL, c, 0},
    {"d", NULL, d, 0},
  };

  gsl_function model_f;
  model_f.function = &lotka_volterra_wrap;
  model_f.params = &params;
  (void)model_f;

  gsl_function kernel_f;
  kernel_f.function = &kde;

  gsl_matrix* params_weighted = NULL;
  params_weighted = gsl_matrix_alloc(2, 4);

  gsl_function sampling;
  gsl_function proposal;

  int ret = abc_smc(params_weighted,
		    &model_f,
		    &kernel_f,
		    2,
		    &sampling,
		    &proposal,
		    1);
  if (ret != GSL_SUCCESS) {
    GSL_ERROR("Sample failed", ret);
  }

  gsl_matrix_free(params_weighted);

  return EXIT_SUCCESS;
}

