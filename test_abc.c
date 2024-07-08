/**
   @file test_abc.c

   @brief Unit test suite for MPI-ABC public API.
 */

#include <abc.h>

#include <check.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>


/**
   Unit test ODE model.
 */
START_TEST(test_ode)
{
  /* Allocate model params vector and outcome matrix. */
  gsl_vector* params = NULL;
  params = gsl_vector_alloc(4);
  ck_assert_ptr_nonnull(params);
  gsl_vector_set(params, 0, 1.00); /* a. */
  gsl_vector_set(params, 1, 1.00); /* b (NB: detuned from 0.1). */
  gsl_vector_set(params, 2, 1.50); /* c. */
  gsl_vector_set(params, 3, 0.75); /* d. */
  const size_t T = 100;
  const size_t n_eq = 2;
  gsl_matrix* outcomes = NULL;
  outcomes = gsl_matrix_alloc(T, 1 + n_eq);
  ck_assert_ptr_nonnull(outcomes);

  /* Calculate ODE equation solutions. */
  int ret = lotka_volterra_run(outcomes, params);
  ck_assert_int_eq(ret, GSL_SUCCESS);
  double v0 = lotka_volterra_verify(10.0, 5.0, params);
  /* Verify the solution is valid using the V implicit relationship. */
  for (size_t i = 0; i < T; ++i) {
    //double t = gsl_matrix_get(outcomes, i, 0);
    double y0 = gsl_matrix_get(outcomes, i, 1);
    double y1 = gsl_matrix_get(outcomes, i, 2);
    double v = lotka_volterra_verify(y0, y1, params);
    //printf("%.5e %.5e %.5e %.5e\n", t, y0, y1, v);
    ck_assert_double_eq_tol(v, v0, 1e-5);
  }

  /* Cleanup. */
  gsl_vector_free(params);
}


/**
   Unit test the kernel density estimate.
 */
START_TEST(test_kde)
{
  /* Populate a bimodal gaussian distribution using the Wikipedia example:
     https://en.wikipedia.org/wiki/Kernel_density_estimation#Example */
  gsl_vector* data = NULL;
  data = gsl_vector_alloc(6);
  ck_assert_ptr_nonnull(data);
  gsl_vector_set(data, 0, -2.1);
  gsl_vector_set(data, 1, -1.3);
  gsl_vector_set(data, 2, -0.4);
  gsl_vector_set(data, 3,  1.9);
  gsl_vector_set(data, 4,  5.1);
  gsl_vector_set(data, 5,  6.2);

  /* Integrate the kde function. */
  double h = silverman(data);
  //printf("Silverman rule of thumb (h): %.5e\n", h);
  ck_assert_double_gt(h, 0);
  kde_params_t kde_params = {data, h};
  gsl_function kde_f;
  kde_f.function = &kde;
  kde_f.params = &kde_params;
  double a = gsl_vector_get(data, 0) - 3 * h;
  double b = gsl_vector_get(data, 5) + 3 * h;
  double sum = 0;
  double abserr = 0;
  size_t neval = 0;
  int ret =
    gsl_integration_qng(&kde_f, a, b, 1e-3, 1e-3, &sum, &abserr, &neval);
  ck_assert_int_eq(ret, GSL_SUCCESS);
  ck_assert_double_eq_tol(sum, 1.0, 0.01);
  //printf("Integrated KDE: sum = %.5e, abserr = %.5e, neval = %zu\n",
  //       sum, abserr, neval);

  /* Cleanup. */
  gsl_vector_free(data);
}
END_TEST


/**
   Aggregate all unit tests into a suite.
 */
Suite*
suite_abc(void)
{
  Suite *s = suite_create("abc");
  TCase *tc_model = tcase_create("model");
  tcase_add_test(tc_model, test_ode);
  suite_add_tcase(s, tc_model);
  TCase *tc_sampler = tcase_create("sampler");
  tcase_add_test(tc_sampler, test_kde);
  suite_add_tcase(s, tc_sampler);
  return s;
}


/**
   Run unit tests.
 */
int
main(void)
{
  int n_failed = 0;
  Suite *s = suite_abc();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  n_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (n_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
