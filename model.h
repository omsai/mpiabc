/*
  Copyright (c) 2024, Pariksheet Nanda

  SPDX-License-Identifier: Apache-2.0
*/
#ifndef MPIABC_MODEL_H
#define MPIABC_MODEL_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>


/** Common function signatures of gsl_ran_*(...). */
typedef union {
  double(*func0p)(const gsl_rng*);
  double(*func1p)(const gsl_rng*, const double);
  double(*func2p)(const gsl_rng*, const double, const double);
} gsl_ran_function;

/** Prior paramater distributions and fixed parameters. */
typedef struct {
  char* name;
  gsl_ran_function* distr;	/* Use NULL for fixed parameters. */
  double* distr_hyperparams;
  size_t n_hyperparams;
} param_t;


/* Public API. */
double param_sample(const gsl_rng *r, param_t* param);
double lotka_volterra_wrap(double x, void* params_);

/* Unit testing API. */
int lotka_volterra_run(gsl_matrix* outcomes, gsl_vector* params);
double lotka_volterra_verify(double y0, double y1, gsl_vector* params);

#endif
