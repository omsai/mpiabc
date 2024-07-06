/**
   @file test_abc.c

   @brief Unit test suite for MPI-ABC public API.
 */

#include <abc.h>

#include <check.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>


/**
   Unit test the kernel density estimate.
 */
START_TEST(test_kde)
{
  /* Populate a bimodal gaussian distribution using the Wikipedia example:
     https://en.wikipedia.org/wiki/Kernel_density_estimation#Example */
  gsl_block* data = NULL;
  data = gsl_block_alloc(6);
  ck_assert_ptr_nonnull(data);
  data->data[0] = -2.1;
  data->data[1] = -1.3;
  data->data[2] = -0.4;
  data->data[3] =  1.9;
  data->data[4] =  5.1;
  data->data[5] =  6.2;

  /* Call the kde function. */
  double h = silverman(data);
  //printf("Silverman rule of thumb (h): %.5e\n", h);
  gsl_block *x = NULL;
  x = gsl_block_alloc(100);
  ck_assert_ptr_nonnull(x);
  double x_begin = data->data[0] - 3 * h;
  double x_end   = data->data[5] + 3 * h;
  for (size_t i = 0; i < x->size; ++i) {
    x->data[i] = x_begin + (x_end - x_begin) * i / (x->size - 1);
  }
  gsl_block *y = NULL;
  int ret = kde(&y, x, data, h);

  /* Check the returned values. */
  ck_assert_int_eq(ret, GSL_SUCCESS);
  ck_assert_ptr_nonnull(y);
  gsl_vector *y_integral = NULL;
  y_integral = gsl_vector_alloc(y->size - 1);
  for (size_t i = 0; i < y->size - 1; ++i) {
    y_integral->data[i] =
      (y->data[i+1] + y->data[i]) / 2 *
      (x->data[i+1] - x->data[i]);
  }
  ck_assert_ptr_nonnull(y_integral);
  double sum = gsl_vector_sum(y_integral);
  //printf("Integrated KDE sum: %.5e\n", sum);
  ck_assert_double_le(sum, 1);
  ck_assert_double_eq_tol(sum, 1.0, 0.01);

  /* Export kde.gnuplot file; visualize with `gnuplot -p kde.gnuplot`. */
  FILE* gnuplot = fopen("kde.gnuplot", "w");
  fprintf(gnuplot, "unset key\n");
  fprintf(gnuplot, "array data[6] = [%+3.1e", data->data[0]);
  for (size_t i = 1; i < data->size; ++i) {
    fprintf(gnuplot, ", %+3.1f", data->data[i]);
  }
  fprintf(gnuplot, "]\narray zeros[6] = [0, 0, 0, 0, 0, 0]\n");
  fprintf(gnuplot, "$kde << EOD\n");
  for (size_t i = 0; i < x->size; ++i) {
    fprintf(gnuplot, "%+.5e %.5e\n", x->data[i], y->data[i]);
  }
  fprintf(gnuplot, "EOD\n");
  fprintf(gnuplot,
	  "plot sample [i=1:6] '+' using (data[i]):(zeros[i]) \\\n"
	  "     with points ps 2 pt 7, \\\n"
	  "     \"$kde\" with lines \\\n"
	  "     title \"kde\"\n");
  fflush(gnuplot);

  /* Cleanup. */
  gsl_vector_free(y_integral);
  gsl_block_free(y);
  gsl_block_free(x);
  gsl_block_free(data);
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
