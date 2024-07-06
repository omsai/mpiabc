/**
   @file test_abc.c

   @brief Unit test suite for MPI-ABC public API.
 */

#include <abc.h>

#include <check.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>


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
   Aggregate sampler unit tests into a suite.
 */
Suite*
suite_abc(void)
{
  Suite *s = suite_create("abc");
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
