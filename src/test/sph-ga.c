#include <stdio.h>
#include <string.h>
#include <sph/status.h>
#include <sph/test.h>
#include <sph-ga/sph-ga.h>

typedef struct {
  char *input;
  char *expected;
} ga_unary_string_case_t;

typedef struct {
  char *left;
  char *right;
  char *expected;
} ga_binary_string_case_t;

typedef struct {
  char *input;
  char *expected;
} ga_roundtrip_string_case_t;

ga_t test_ga_r3;
ga_t test_ga_c3;

status_t test_ga_run_roundtrip_cases(char *name, ga_t *ga, ga_roundtrip_string_case_t *case_array, uint64_t case_count)
{
  status_declare;
  ga_term_t mv_storage[ga_blade_count];
  ga_mv_t mv;
  char output[ga_mv_string_capacity];
  uint64_t case_index;

  mv.term = mv_storage;
  mv.term_count = 0;

  case_index = 0;
  while (case_index < case_count) {
    ga_mv_from_string(ga, case_array[case_index].input, (&mv));
    ga_mv_to_string(ga, (&mv), output);
    test_helper_assert3(
      name,
      (!strcmp(output, case_array[case_index].expected)),
      "%s", case_array[case_index].input,
      "%s", case_array[case_index].expected,
      "%s", output
    );
    case_index += 1;
  }

exit:
  status_return;
}

status_t test_ga_run_unary_cases(char *name, ga_t *ga, void (*operator)(ga_t *ga, ga_mv_t *in, ga_mv_t *out), ga_unary_string_case_t *case_array, uint64_t case_count)
{
  status_declare;
  ga_term_t input_storage[ga_blade_count];
  ga_term_t output_storage[ga_blade_count];
  ga_mv_t input;
  ga_mv_t output;
  char actual[ga_mv_string_capacity];
  uint64_t case_index;

  input.term = input_storage;
  input.term_count = 0;
  output.term = output_storage;
  output.term_count = 0;

  case_index = 0;
  while (case_index < case_count) {
    ga_mv_from_string(ga, case_array[case_index].input, (&input));
    operator(ga, (&input), (&output));
    ga_mv_to_string(ga, (&output), actual);
    test_helper_assert3(
      name,
      (!strcmp(actual, case_array[case_index].expected)),
      "%s", case_array[case_index].input,
      "%s", case_array[case_index].expected,
      "%s", actual
    );
    case_index += 1;
  }

exit:
  status_return;
}

status_t test_ga_run_binary_cases(char *name, ga_t *ga, status_t (*operator)(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out), ga_binary_string_case_t *case_array, uint64_t case_count)
{
  status_declare;
  ga_term_t left_storage[ga_blade_count];
  ga_term_t right_storage[ga_blade_count];
  ga_term_t output_storage[ga_blade_count];
  ga_mv_t left;
  ga_mv_t right;
  ga_mv_t output;
  char actual[ga_mv_string_capacity];
  uint64_t case_index;

  left.term = left_storage;
  left.term_count = 0;
  right.term = right_storage;
  right.term_count = 0;
  output.term = output_storage;
  output.term_count = 0;

  case_index = 0;
  while (case_index < case_count) {
    ga_mv_from_string(ga, case_array[case_index].left, (&left));
    ga_mv_from_string(ga, case_array[case_index].right, (&right));
    status_require((operator(ga, (&left), (&right), (&output))));
    ga_mv_to_string(ga, (&output), actual);
    test_helper_assert4(
      name,
      (!strcmp(actual, case_array[case_index].expected)),
      "%s", case_array[case_index].left,
      "%s", case_array[case_index].right,
      "%s", case_array[case_index].expected,
      "%s", actual
    );
    case_index += 1;
  }

exit:
  status_return;
}

status_t test_ga_string_roundtrip_c3(void)
{
  status_declare;
  ga_roundtrip_string_case_t case_array[] = {
    {"0", "0"},
    {"1", "1"},
    {"-1", "-1"},
    {"e1", "e1"},
    {"-e1", "-e1"},
    {"e1_2", "e1_2"},
    {"2e1_2", "2e1_2"},
    {"1 + e1", "1 + e1"},
    {"1 + e1 + -e1_2", "1 + e1 + -e1_2"},
    {"54e1_2_3_5 + 48e1_4_5 + 63e2_3_4_5", "54e1_2_3_5 + 48e1_4_5 + 63e2_3_4_5"},
    {"no", "e4"},
    {"ni", "e5"},
    {"1 + no + -2ni", "1 + e4 + -2e5"}
  };

  status_require((test_ga_run_roundtrip_cases("ga string roundtrip c3", (&test_ga_c3), case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_reverse(void)
{
  status_declare;
  ga_unary_string_case_t case_array[] = {
    {"e1", "e1"},
    {"e1_2", "-e1_2"},
    {"e1_2_3", "-e1_2_3"},
    {"1 + e1", "1 + e1"},
    {"1 + e1 + e1_2", "1 + e1 + -e1_2"},
    {"1 + e1 + e2 + e1_2 + e1_3 + e1_2_3", "1 + e1 + e2 + -e1_2 + -e1_3 + -e1_2_3"}
  };

  status_require((test_ga_run_unary_cases("ga reverse", (&test_ga_c3), ga_reverse, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_involute(void)
{
  status_declare;
  ga_unary_string_case_t case_array[] = {
    {"e1", "-e1"},
    {"e1_2", "e1_2"},
    {"e1_2_3", "-e1_2_3"}
  };

  status_require((test_ga_run_unary_cases("ga involute", (&test_ga_c3), ga_involute, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_conjugate(void)
{
  status_declare;
  ga_unary_string_case_t case_array[] = {
    {"e1", "-e1"},
    {"e1_2", "-e1_2"},
    {"e1_2_3", "e1_2_3"}
  };

  status_require((test_ga_run_unary_cases("ga conjugate", (&test_ga_c3), ga_conjugate, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_gp_c3(void)
{
  status_declare;
  ga_binary_string_case_t case_array[] = {
    {"e1_2_3_4_5", "e1_2_3_4_5", "1"},
    {"2", "3", "6"},
    {"7", "8e1", "56e1"},
    {"4", "5e1_2", "20e1_2"},
    {"5e2", "6", "30e2"},
    {"2e3", "3e3", "6"},
    {"9e1", "8e2", "72e1_2"},
    {"7e1_3", "8e3_5", "56e1_5"},
    {"9e1_2", "8e1_2", "-72"},
    {"4e1", "5e5", "20e1_5"},
    {"7e3", "8e1_3", "-56e1"},
    {"9e1", "2e1_5", "18e5"},
    {"4e3", "5e1_2_3", "20e1_2"},
    {"6e1", "7e1_3_4", "42e3_4"},
    {"8e4", "9", "72e4"},
    {"5e4", "6e2", "-30e2_4"},
    {"7e5", "8e5", "0"},
    {"2e5", "3e1_3", "6e1_3_5"},
    {"4e4", "5e2_4", "0"},
    {"7e4", "8e1_2_3", "-56e1_2_3_4"},
    {"3e2_3", "4", "12e2_3"},
    {"9e5", "2e2_4_5", "-18e2_5"},
    {"8e1_2", "9e1", "-72e2"},
    {"2e1_3", "3e4", "6e1_3_4"},
    {"5e1_2", "6e1_3", "-30e2_3"},
    {"2e1_2", "3e1_2_3", "-6e3"},
    {"6e2_4", "7", "42e2_4"},
    {"4e1_3", "5e1_3_5", "-20e5"},
    {"3e3_5", "4e3", "-12e5"},
    {"5e4_5", "6e5", "-30e5"},
    {"7e3_5", "8e1_3_5", "0"},
    {"2e3_4", "3e4_5", "-6e3_4"},
    {"8e2_4", "9e2_3", "72e3_4"},
    {"5e2_5", "6e1_2_3", "30e1_3_5"},
    {"6e4 + 7e5", "8e2_5 + 9e2_3_5", "48e2 + -54e2_3"}
  };

  status_require((test_ga_run_binary_cases("ga gp c3", (&test_ga_c3), ga_gp_one, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_gp_r3(void)
{
  status_declare;
  ga_binary_string_case_t case_array[] = {
    {"e1_2_3", "e1_2_3", "-1"}
  };

  status_require((test_ga_run_binary_cases("ga gp r3", (&test_ga_r3), ga_gp_one, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_ip_c3(void)
{
  status_declare;
  ga_binary_string_case_t case_array[] = {
    {"e2_4", "e2_4_5", "e4"},
    {"e2_4_5", "e2_4_5", "-1"},
    {"e4 + e1_2 + e3_5", "e1 + e2_3 + e4_5", "e4"},
    {"e4", "e5", "-1"},
    {"e5", "e4", "-1"},
    {"e4 + e5", "e4 + e5", "-2"},
    {"6e1 + 7e4", "8e4_5 + 9e2_3_5", "-63e2_3 + 56e4"},
    {"4e4_1_2", "3e4_1_2_5", "12e4"},
    {"e1_4", "e1_5", "-1"},
    {"e1_5", "e1_4", "-1"},
    {"2e1_2_3", "e1_2_3_4_5", "2e4_5"},
    {"4", "5e1", "20e1"},
    {"6e2", "7", "0"},
    {"8e1_2", "9e1", "0"},
    {"9e5", "2e2_4_5", "18e2_5"},
    {"e5", "e5", "0"},
    {"e4", "e4", "0"},
    {"4e1", "5e1", "20"},
    {"3e4", "4e1_3_4", "0"},
    {"8e3", "9e1", "0"},
    {"2e2", "3e2", "6"},
    {"4e3", "5e1_2", "0"},
    {"7e1_3", "8e3_5", "0"},
    {"4e2", "5e2_3", "20e3"},
    {"7e2", "8e4", "0"},
    {"8e2", "9e1_2_4", "-72e1_4"},
    {"3e5", "4e1", "0"},
    {"6e5", "7e1_3", "0"},
    {"9e5", "2e4", "-18"},
    {"3e5", "4e5", "0"},
    {"5e4", "6e4_5", "30e4"},
    {"e4", "e4", "0"},
    {"e5", "e5", "0"},
    {"e4", "e5", "-1"},
    {"e5", "e4", "-1"},
    {"e4", "e1", "0"},
    {"e1", "e5", "0"},
    {"e2", "e3", "0"},
    {"e3", "e3", "1"},
    {"e4_5", "e4_5", "-1"},
    {"e1_4", "e1_4", "0"},
    {"e1_5", "e1_5", "0"},
    {"e4_5", "e1_4", "0"},
    {"e4_5", "e1_5", "0"},
    {"e4_5", "e2_3", "0"},
    {"e1_4_5", "e1_4_5", "-1"},
    {"e1_2_5", "e1_4_5", "0"},
    {"e1_2_4_5", "e1_2_4_5", "1"},
    {"e1_2_4_5", "e1_2_3_5", "0"},
    {"e1_2_3_4_5", "e1_2_3_4_5", "1"},
    {"e1 + e5", "e1 + e5", "1"},
    {"2", "3e4", "6e4"},
    {"5", "7e1_2", "35e1_2"}
  };

  status_require((test_ga_run_binary_cases("ga ip c3", (&test_ga_c3), ga_ip_one, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_ep_c3(void)
{
  status_declare;
  ga_binary_string_case_t case_array[] = {
    {"2", "3", "0"},
    {"4", "5e1", "20e1"},
    {"6e2", "7", "42e2"},
    {"3e4", "4e1_3_4", "0"},
    {"7e5", "8e5", "0"},
    {"7e1_3", "8e3_5", "0"},
    {"8e3", "9e1", "-72e1_3"},
    {"2e2", "3e2", "0"},
    {"4e3", "5e1_2", "20e1_2_3"},
    {"4e2", "5e2_3", "0"},
    {"7e2", "8e4", "56e2_4"},
    {"2e1", "3e4_5", "6e1_4_5"},
    {"6e1_4", "5e3", "-30e1_3_4"},
    {"8e2", "9e1_2_4", "0"},
    {"9e1_5", "2e3", "-18e1_3_5"},
    {"6e1 + 7e4", "8e4_5 + 9e2_3_5", "54e1_2_3_5 + 48e1_4_5 + 63e2_3_4_5"}
  };

  status_require((test_ga_run_binary_cases("ga ep c3", (&test_ga_c3), ga_ep_one, case_array, ((uint64_t) (sizeof(case_array) / sizeof(case_array[0]))))));

exit:
  status_return;
}

status_t test_ga_no_and_ni(void)
{
  status_declare;
  ga_term_t out_storage[ga_blade_count];
  ga_mv_t out;
  char actual[ga_mv_string_capacity];

  out.term = out_storage;
  out.term_count = 0;

  ga_no((&test_ga_c3), 2.0, (&out));
  ga_mv_to_string((&test_ga_c3), (&out), actual);
  test_helper_assert1("ga no", (!strcmp(actual, "2e4")), "%s", actual);

  ga_ni((&test_ga_c3), -3.0, (&out));
  ga_mv_to_string((&test_ga_c3), (&out), actual);
  test_helper_assert1("ga ni", (!strcmp(actual, "-3e5")), "%s", actual);

exit:
  status_return;
}

status_t test_ga_point_and_point_euclidean(void)
{
  status_declare;
  ga_term_t point_storage[ga_blade_count];
  ga_mv_t point;
  char actual[ga_mv_string_capacity];
  double euclidean_input[3];
  double euclidean_output[3];

  point.term = point_storage;
  point.term_count = 0;

  euclidean_input[0] = 2.0;
  euclidean_input[1] = 3.0;
  euclidean_input[2] = 4.0;

  ga_point((&test_ga_c3), euclidean_input, (&point));
  ga_mv_to_string((&test_ga_c3), (&point), actual);
  test_helper_assert1("ga point string", (!strcmp(actual, "2e1 + 3e2 + 4e3 + e4 + 14.5e5")), "%s", actual);

  ga_point_euclidean((&test_ga_c3), (&point), euclidean_output);
  test_helper_assert1("ga point euclidean x", (euclidean_output[0] == 2.0), "%f", euclidean_output[0]);
  test_helper_assert1("ga point euclidean y", (euclidean_output[1] == 3.0), "%f", euclidean_output[1]);
  test_helper_assert1("ga point euclidean z", (euclidean_output[2] == 4.0), "%f", euclidean_output[2]);

exit:
  status_return;
}

status_t test_ga_rotor(void)
{
  status_declare;
  ga_term_t rotor_storage[ga_blade_count];
  ga_mv_t rotor;
  char actual[ga_mv_string_capacity];
  double coeff[4];

  rotor.term = rotor_storage;
  rotor.term_count = 0;

  coeff[0] = 2.0;
  coeff[1] = 3.0;
  coeff[2] = 4.0;
  coeff[3] = 5.0;

  ga_rotor((&test_ga_c3), coeff, (&rotor));
  ga_mv_to_string((&test_ga_c3), (&rotor), actual);
  test_helper_assert1("ga rotor string", (!strcmp(actual, "2 + 3e1_2 + 4e1_3 + 5e2_3")), "%s", actual);

exit:
  status_return;
}

int main(void)
{
  status_declare;
  double r3_metric[3];
  double c3_metric[3];
  uint64_t r3_inited;
  uint64_t c3_inited;

  r3_metric[0] = 1.0;
  r3_metric[1] = 1.0;
  r3_metric[2] = 1.0;

  c3_metric[0] = 1.0;
  c3_metric[1] = 1.0;
  c3_metric[2] = 1.0;

  r3_inited = 0;
  c3_inited = 0;

  test_ga_r3.dimensions = 0;
  test_ga_c3.dimensions = 0;

  status_require((ga_init_diagonal((&test_ga_r3), r3_metric, 3, 0)));
  r3_inited = 1;
  status_require((ga_init_diagonal((&test_ga_c3), c3_metric, 3, 1)));
  c3_inited = 1;

  test_helper_test(test_ga_string_roundtrip_c3);
  test_helper_test(test_ga_reverse);
  test_helper_test(test_ga_involute);
  test_helper_test(test_ga_conjugate);
  test_helper_test(test_ga_gp_c3);
  test_helper_test(test_ga_gp_r3);
  test_helper_test(test_ga_ip_c3);
  test_helper_test(test_ga_ep_c3);
  test_helper_test(test_ga_no_and_ni);
  test_helper_test(test_ga_point_and_point_euclidean);
  test_helper_test(test_ga_rotor);

exit:
  if (c3_inited) {
    ga_uninit((&test_ga_c3));
  }
  if (r3_inited) {
    ga_uninit((&test_ga_r3));
  }
  test_helper_display_summary_description(ga_status_description);
  return status.id;
}