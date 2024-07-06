#ifndef MPIABC_SAMPLER_H
#define MPIABC_SAMPLER_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>


typedef struct {
  gsl_vector* data;
  double bandwidth;
} kde_params_t;


double silverman(gsl_vector* x);
double kde(double x, void* params_);

#endif
