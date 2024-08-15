/*
  Copyright (c) 2024, Pariksheet Nanda

  SPDX-License-Identifier: Apache-2.0
*/
#ifndef MPIABC_SAMPLER_H
#define MPIABC_SAMPLER_H

#include <model.h>


/** Convenience to run gsl_integration_qng() on kde(). */
typedef struct {
  gsl_vector* data;
  double bandwidth;
} kde_params_t;


double silverman(gsl_vector* x);
double kde(double x, void* params_);

int abc_smc(gsl_matrix* params,
	    gsl_function* model,
	    gsl_function* kernel,
	    int n,
	    gsl_function* sampling,
	    gsl_function* proposal,
	    int alpha);

#endif
