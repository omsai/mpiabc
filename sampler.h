#ifndef SAMPLER_H
#define SAMPLER_H

#include <gsl/gsl_block.h>

double silverman(gsl_block* x);
int kde(gsl_block** y, gsl_block* x, gsl_block* data, double h);

#endif
