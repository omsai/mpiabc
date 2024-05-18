#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_odeiv2.h>


int
lotka_volterra(double t, const double y[], double dydt[], void* params_) {
    (void)(t);			/* Suppress unused parameter warning. */
    gsl_block* params = params_;
    double a = params->data[0];
    double b = params->data[1];
    double c = params->data[2];
    double d = params->data[3];

    dydt[0] =  a * y[0] -     b * y[0] * y[1];
    dydt[1] = -c * y[1] + d * b * y[0] * y[1];

    return GSL_SUCCESS;
}


int
main(void) {
    gsl_block* params = gsl_block_alloc(4);
    if (params == NULL) {
	GSL_ERROR("could not allocate params", ENOMEM);
    }
    params->data[0] = 1.0;	/* a. */
    params->data[1] = 1.0;	/* b (NB: detuned from 0.1). */
    params->data[2] = 1.5;	/* c. */
    params->data[3] = 0.75;	/* d. */
    gsl_odeiv2_system model = {
	lotka_volterra,		/* ODE function. */
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
    
    for (int i = 1; i <= 100; ++i) {
	double ti = i * t_end / 100.0;
	int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
	if (status != GSL_SUCCESS) {
	    gsl_odeiv2_driver_free(driver);
	    gsl_block_free(params);
	    GSL_ERROR("ODE failed", status);
	}
	printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

    gsl_odeiv2_driver_free(driver);
    gsl_block_free(params);
    return 0;
}
