#include <stdio.h>
#include <sph/status.h>
#include <sph/test.h>
#include <sph-ga/sph-ga.h>

status_t test_ga_init_diagonal(void)
{
  status_declare;
  ga_t ga;
  double metric[3];

  metric[0] = 1.0;
  metric[1] = 1.0;
  metric[2] = 1.0;
  ga.dimensions = 0;

  status_require((ga_init_diagonal((&ga), metric, 3, 0)));
  test_helper_assert1("ga init diagonal dimensions", (ga.dimensions == 3), "%" PRId64, ga.dimensions);
  test_helper_assert("ga init diagonal symmetric", (ga.is_symmetric));
  test_helper_assert("ga init diagonal diagonal", (ga.is_diagonal));
  test_helper_assert1("ga init diagonal pseudoscalar id", (ga.pseudoscalar_id == 7), "%" PRIu64, ((uint64_t) ga.pseudoscalar_id));

exit:
  if (ga.dimensions) {
    ga_uninit((&ga));
  }
  status_return;
}

status_t test_ga_ep_basis_antisymmetry(void)
{
  status_declare;
  ga_t ga;
  double metric[2];
  ga_term_t left_storage[1];
  ga_term_t right_storage[1];
  ga_term_t out_a_storage[ga_blade_count];
  ga_term_t out_b_storage[ga_blade_count];
  ga_mv_t left;
  ga_mv_t right;
  ga_mv_t out_a;
  ga_mv_t out_b;
  ga_accumulator_t accumulator_a;
  ga_accumulator_t accumulator_b;

  metric[0] = 1.0;
  metric[1] = 1.0;
  ga.dimensions = 0;
  accumulator_a.used_count = 0;
  accumulator_b.used_count = 0;

  left.term = left_storage;
  right.term = right_storage;
  out_a.term = out_a_storage;
  out_b.term = out_b_storage;

  status_require((ga_init_diagonal((&ga), metric, 2, 0)));

  ga_basis_blade(1, 1.0, (&left));
  ga_basis_blade(2, 1.0, (&right));

  ga_ep_one((&ga), (&left), (&right), (&accumulator_a), (&out_a));
  ga_ep_one((&ga), (&right), (&left), (&accumulator_b), (&out_b));

  test_helper_assert1("ga ep basis term count a", (out_a.term_count == 1), "%" PRId64, out_a.term_count);
  test_helper_assert1("ga ep basis term count b", (out_b.term_count == 1), "%" PRId64, out_b.term_count);
  test_helper_assert1("ga ep basis id a", (out_a.term[0].id == 3), "%" PRIu64, ((uint64_t) out_a.term[0].id));
  test_helper_assert1("ga ep basis id b", (out_b.term[0].id == 3), "%" PRIu64, ((uint64_t) out_b.term[0].id));
  test_helper_assert1("ga ep basis coeff a", (out_a.term[0].coeff == 1.0), "%f", out_a.term[0].coeff);
  test_helper_assert1("ga ep basis coeff b", (out_b.term[0].coeff == -1.0), "%f", out_b.term[0].coeff);

exit:
  if (ga.dimensions) {
    ga_uninit((&ga));
  }
  status_return;
}

status_t test_ga_gp_basis_square(void)
{
  status_declare;
  ga_t ga;
  double metric[2];
  ga_term_t vector_storage[1];
  ga_term_t out_storage[ga_blade_count];
  ga_mv_t vector;
  ga_mv_t out;
  ga_accumulator_t accumulator;

  metric[0] = 1.0;
  metric[1] = -1.0;
  ga.dimensions = 0;
  accumulator.used_count = 0;

  vector.term = vector_storage;
  out.term = out_storage;

  status_require((ga_init_diagonal((&ga), metric, 2, 0)));

  ga_basis_blade(1, 1.0, (&vector));
  ga_gp_one((&ga), (&vector), (&vector), (&accumulator), (&out));

  test_helper_assert1("ga gp e1 square term count", (out.term_count == 1), "%" PRId64, out.term_count);
  test_helper_assert1("ga gp e1 square id", (!out.term[0].id), "%" PRIu64, ((uint64_t) out.term[0].id));
  test_helper_assert1("ga gp e1 square coeff", (out.term[0].coeff == 1.0), "%f", out.term[0].coeff);

  ga_basis_blade(2, 1.0, (&vector));
  ga_gp_one((&ga), (&vector), (&vector), (&accumulator), (&out));

  test_helper_assert1("ga gp e2 square term count", (out.term_count == 1), "%" PRId64, out.term_count);
  test_helper_assert1("ga gp e2 square id", (!out.term[0].id), "%" PRIu64, ((uint64_t) out.term[0].id));
  test_helper_assert1("ga gp e2 square coeff", (out.term[0].coeff == -1.0), "%f", out.term[0].coeff);

exit:
  if (ga.dimensions) {
    ga_uninit((&ga));
  }
  status_return;
}

status_t test_ga_inverse_vector(void)
{
  status_declare;
  ga_t ga;
  double metric[2];
  ga_term_t vector_storage[1];
  ga_term_t reverse_storage[ga_blade_count];
  ga_term_t denominator_storage[ga_blade_count];
  ga_term_t inverse_storage[ga_blade_count];
  ga_term_t product_storage[ga_blade_count];
  ga_mv_t vector;
  ga_mv_t reverse_mv;
  ga_mv_t denominator_mv;
  ga_mv_t inverse_mv;
  ga_mv_t product_mv;
  ga_accumulator_t accumulator_inverse;
  ga_accumulator_t accumulator_product;

  metric[0] = 1.0;
  metric[1] = 1.0;
  ga.dimensions = 0;
  accumulator_inverse.used_count = 0;
  accumulator_product.used_count = 0;

  vector.term = vector_storage;
  reverse_mv.term = reverse_storage;
  denominator_mv.term = denominator_storage;
  inverse_mv.term = inverse_storage;
  product_mv.term = product_storage;

  status_require((ga_init_diagonal((&ga), metric, 2, 0)));

  ga_basis_blade(1, 2.0, (&vector));
  status_require((ga_inverse((&ga), (&vector), (&accumulator_inverse), (&reverse_mv), (&denominator_mv), (&inverse_mv))));
  ga_gp_one((&ga), (&vector), (&inverse_mv), (&accumulator_product), (&product_mv));

  test_helper_assert1("ga inverse vector term count", (inverse_mv.term_count == 1), "%" PRId64, inverse_mv.term_count);
  test_helper_assert1("ga inverse vector id", (inverse_mv.term[0].id == 1), "%" PRIu64, ((uint64_t) inverse_mv.term[0].id));
  test_helper_assert1("ga inverse vector coeff", (inverse_mv.term[0].coeff == 0.5), "%f", inverse_mv.term[0].coeff);

  test_helper_assert1("ga inverse product term count", (product_mv.term_count == 1), "%" PRId64, product_mv.term_count);
  test_helper_assert1("ga inverse product id", (!product_mv.term[0].id), "%" PRIu64, ((uint64_t) product_mv.term[0].id));
  test_helper_assert1("ga inverse product coeff", (product_mv.term[0].coeff == 1.0), "%f", product_mv.term[0].coeff);

exit:
  if (ga.dimensions) {
    ga_uninit((&ga));
  }
  status_return;
}

status_t test_ga_init_conformal(void)
{
  status_declare;
  ga_t ga;
  double metric[3];

  metric[0] = 1.0;
  metric[1] = 1.0;
  metric[2] = 1.0;
  ga.dimensions = 0;

  status_require((ga_init_diagonal((&ga), metric, 3, 1)));

  test_helper_assert1("ga conformal dimensions", (ga.dimensions == 5), "%" PRId64, ga.dimensions);
  test_helper_assert("ga conformal symmetric", (ga.is_symmetric));
  test_helper_assert("ga conformal diagonal classification", (!ga.is_diagonal));
  test_helper_assert1("ga conformal null count", (ga.null_vector_count == 2), "%" PRId64, ga.null_vector_count);
  test_helper_assert1("ga conformal no id", (ga.no_id == 8), "%" PRIu64, ((uint64_t) ga.no_id));
  test_helper_assert1("ga conformal ni id", (ga.ni_id == 16), "%" PRIu64, ((uint64_t) ga.ni_id));

exit:
  if (ga.dimensions) {
    ga_uninit((&ga));
  }
  status_return;
}

int main(void)
{
  status_declare;

  test_helper_test(test_ga_init_diagonal);
  test_helper_test(test_ga_ep_basis_antisymmetry);
  test_helper_test(test_ga_gp_basis_square);
  test_helper_test(test_ga_inverse_vector);
  test_helper_test(test_ga_init_conformal);

exit:
  test_helper_display_summary_description(ga_status_description);
  return status.id;
}
