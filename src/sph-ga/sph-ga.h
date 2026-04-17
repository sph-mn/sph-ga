#ifndef sph_ga_h_included
#define sph_ga_h_included

#include <inttypes.h>
#include <sph/status.h>

#ifndef ga_max_dimensions
#define ga_max_dimensions 10
#endif

#ifndef ga_blade_id_t
#define ga_blade_id_t uint64_t
#endif

#define ga_blade_count ((((ga_blade_id_t) 1)) << ga_max_dimensions)

enum {
  ga_status_id_metric_not_symmetric = 1,
  ga_status_id_conformal_requires_two_null_vectors = 2,
  ga_status_id_conformal_requires_diagonal_metric = 3,
  ga_status_id_inverse_zero_denominator = 4
};

typedef struct {
  ga_blade_id_t id;
  double coeff;
} ga_term_t;

typedef struct {
  ga_term_t *term;
  int64_t term_count;
} ga_mv_t;

typedef struct {
  double coeff[ga_blade_count];
  ga_blade_id_t used_id[ga_blade_count];
  int64_t used_count;
} ga_accumulator_t;

typedef struct {
  int64_t dimensions;
  int64_t is_conformal;
  int64_t is_symmetric;
  int64_t is_diagonal;
  int64_t null_vector_count;
  int64_t null_vector_start;
  ga_blade_id_t id_null;
  ga_blade_id_t pseudoscalar_id;
  ga_blade_id_t no_id;
  ga_blade_id_t ni_id;
  int64_t no_bit_index;
  int64_t ni_bit_index;
  int64_t no_index;
  int64_t ni_index;
  double metric[ga_max_dimensions][ga_max_dimensions];
  uint8_t grade_by_id[ga_blade_count];
} ga_t;

char *ga_status_description(status_t status);

status_t ga_init_full(ga_t *ga, double metric[ga_max_dimensions][ga_max_dimensions], int64_t dimensions, int64_t is_conformal);
status_t ga_init_diagonal(ga_t *ga, double *diagonal_metric, int64_t dimensions, int64_t is_conformal);
void ga_uninit(ga_t *ga);

void ga_accumulator_clear(ga_accumulator_t *accumulator);
void ga_accumulator_add(ga_accumulator_t *accumulator, ga_blade_id_t id, double coeff);
void ga_accumulator_to_mv(ga_t *ga, ga_accumulator_t *accumulator, ga_mv_t *mv);

void ga_zero(ga_mv_t *mv);
void ga_scalar(double coeff, ga_mv_t *mv);
void ga_basis_blade(int64_t index, double coeff, ga_mv_t *mv);

void ga_scale(ga_mv_t *in, double scalar, ga_mv_t *out);
void ga_combine(ga_t *ga, ga_mv_t *left, ga_mv_t *right, double right_scale, ga_accumulator_t *accumulator, ga_mv_t *out);
void ga_add(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out);
void ga_subtract(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out);

void ga_involute(ga_t *ga, ga_mv_t *in, ga_mv_t *out);
void ga_reverse(ga_t *ga, ga_mv_t *in, ga_mv_t *out);
void ga_conjugate(ga_t *ga, ga_mv_t *in, ga_mv_t *out);

void ga_ep_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out);
void ga_gp_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out);
void ga_ip_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out);
status_t ga_inverse(ga_t *ga, ga_mv_t *in, ga_accumulator_t *accumulator, ga_mv_t *reverse_mv, ga_mv_t *denominator_mv, ga_mv_t *out);

#endif
