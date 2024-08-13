/**
   @file model.c

   @brief Model with parameters to calibrate using MPI-ABC.
 */

#include <math.h>
#include <stdio.h>

#include <model.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>


/**
   The classic predator-prey ordinary differential equation (ODE) model.

   The formula used is:
   \f{eqnarray*}{
   {\displaystyle \frac{dx}{dt}} &=& \phantom{-}\alpha x -
                                     \phantom{\delta}\beta xy \\
   {\displaystyle \frac{dy}{dt}} &=& -\gamma y +
                                      \delta \beta xy \\\\
   \mathrm{where,}& \\
   x &\equiv& \textrm{prey population density} \\
   y &\equiv& \textrm{predator population density} \\
   \alpha &\equiv& \textrm{prey natural growth rate} \\
   \beta &\equiv& \textrm{prey interaction rate with predators} \\
   \gamma &\equiv& \textrm{predator natural death rate} \\
   \delta &\equiv& \textrm{predator interaction rate with prey} \\
   \f}

   @param t Time variable necessary for the ODE solver but otherwise unused.
   @param y Equation RHS variables that are treated as read-only.
   @param dydt Equation LHS variables updated by the ODE solver.
   @param params_ Equation parameters stored in a
   [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html).
   @return GSL_SUCCESS
 */
int
lotka_volterra_eqs(double t, const double y[], double dydt[], void* params_) {
    (void)(t);			/* Suppress unused parameter warning. */
    gsl_vector* params = params_;
    double a = gsl_vector_get(params, 0);
    double b = gsl_vector_get(params, 1);
    double c = gsl_vector_get(params, 2);
    double d = gsl_vector_get(params, 3);

    dydt[0] =  a * y[0] -     b * y[0] * y[1];
    dydt[1] = -c * y[1] + d * b * y[0] * y[1];

    return GSL_SUCCESS;
}


/**
   Solve the two predator-prey ordinary differential equations (ODEs).

   @param outcomes Matrix with columns containing time and equation outcomes.
   @param params Equation parameters.

   @return GSL_SUCCESS or gsl_odeiv2_driver_apply() error

   @see lotka_volterra_eqs()
 */
int
lotka_volterra_run(gsl_matrix* outcomes, gsl_vector* params)
{
    gsl_odeiv2_system model = {
	lotka_volterra_eqs,	/* ODE function. */
	NULL,			/* Derivative elements. */
	2,			/* Number of equations. */
	params			/* Pointer to parameters. */
    };
    gsl_odeiv2_driver* driver =
	gsl_odeiv2_driver_alloc_y_new(&model,
				      gsl_odeiv2_step_rkf45,
				      1e-3,  /* Initial step size. */
				      1e-8,  /* epsabs */
				      1e-8); /* epsrel */
    double t = 0.0, t_end = 15.0;
    double y[2] = { 10.0, 5.0 };
    
    for (size_t i = 1; i <= outcomes->size1; ++i) {
	double ti = i * t_end / outcomes->size1;
	int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
	if (status != GSL_SUCCESS) { // GCOVR_EXCL_START
	    gsl_odeiv2_driver_free(driver);
	    GSL_ERROR("ODE failed", status);
	} // GCOVR_EXCL_STOP
	gsl_matrix_set(outcomes, i - 1, 0, t);
	gsl_matrix_set(outcomes, i - 1, 1, y[0]);
	gsl_matrix_set(outcomes, i - 1, 2, y[1]);
    }

    gsl_odeiv2_driver_free(driver);
    return GSL_SUCCESS;
}


/**
   Closed form constant of the implicit relationship between the variables.

   The formula used is:
   \f{eqnarray*}{
   V&=& \delta x - \gamma \ln(x)+ \\
     && \beta y - \alpha \ln(y) \\
   \mathrm{where,}& \\
   x &\equiv& \textrm{prey population density} \\
   y &\equiv& \textrm{predator population density} \\
   \alpha &\equiv& \textrm{prey natural growth rate} \\
   \beta &\equiv& \textrm{prey interaction rate with predators} \\
   \gamma &\equiv& \textrm{predator natural death rate} \\
   \delta &\equiv& \textrm{predator interaction rate with prey} \\
   \f}

   @param y0 Value of the prey variable at a timepoint.
   @param y1 Value of the predator variable at the same timepoint.
   @param params Unchaging parameters for all timepoints.

   @return V constant

   @see lotka_volterra_eqs()
 */
double
lotka_volterra_verify(double y0, double y1, gsl_vector* params)
{
  double a = gsl_vector_get(params, 0);
  double b = gsl_vector_get(params, 1);
  double c = gsl_vector_get(params, 2);
  double d = gsl_vector_get(params, 3);

  return d * y0 - c * log(y0) + b * y1 - a * log(y1);
}


/**
   Return single distance for each parameter.

   @see lotka_volterra_run()
 */
void
lotka_volterra_sum_stat(gsl_matrix* distances, gsl_matrix* outcomes) {
  (void)distances;
  (void)outcomes;
}


/**
   Run model and return distance to the abc_smc() sampler.

   @see abc_smc()
 */
double
lotka_volterra_wrap(double x, void* params_) {
  (void)x;
  (void)params_;
  return EXIT_FAILURE;
}


/**
   Draw a parameter from its probability distribution or its fixed value.

   @param r GSL random number generator to sample from a probability
          distribution.
   @param param Container for fixed vaues or a probability distribution.

   @return A sample from the probability distribution support or NAN if
           param.distr != NULL and param.n_hypermarameters is not within [0, 2].
 */
double
param_sample(const gsl_rng* r, param_t* param) {
  if (param->distr == NULL) {
    return param->distr_hyperparams[0];
  }
  gsl_ran_function* func = param->distr;
  switch (param->n_hyperparams) {
  case 0:
    return func->func0p(r);
  case 1:
    return func->func1p(r,
			param->distr_hyperparams[0]);
  case 2:
    return func->func2p(r,
			param->distr_hyperparams[0],
			param->distr_hyperparams[1]);
  }
  return NAN;
}
