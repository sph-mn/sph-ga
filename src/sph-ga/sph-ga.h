#ifndef sph_ga_h_included
#define sph_ga_h_included

#include <inttypes.h>
#include <sph/status.h>

#ifndef ga_max_dimensions
#define ga_max_dimensions 10
#endif

#ifndef ga_accumulator_stack_size
#define ga_accumulator_stack_size 8
#endif

#ifndef ga_blade_id_t
#define ga_blade_id_t uint64_t
#endif

#ifndef ga_mv_string_capacity
#define ga_mv_string_capacity 512
#endif

#define ga_blade_count ((((ga_blade_id_t) 1)) << ga_max_dimensions)

#define ga_mv_set_zero(mv) \
  do { \
    (mv)->term_count = 0; \
  } while (0)

#define ga_mv_set_scalar(mv, value) \
  do { \
    if ((value)) { \
      (mv)->term[0].id = 0; \
      (mv)->term[0].coeff = (value); \
      (mv)->term_count = 1; \
    } else { \
      (mv)->term_count = 0; \
    } \
  } while (0)

#define ga_mv_set_basis_blade(mv, index, value) \
  do { \
    if (!(value)) { \
      (mv)->term_count = 0; \
    } else if (!(index)) { \
      (mv)->term[0].id = 0; \
      (mv)->term[0].coeff = (value); \
      (mv)->term_count = 1; \
    } else { \
      (mv)->term[0].id = (((ga_blade_id_t) 1)) << ((index) - 1); \
      (mv)->term[0].coeff = (value); \
      (mv)->term_count = 1; \
    } \
  } while (0)

enum {
  ga_status_id_metric_not_symmetric = 1,
  ga_status_id_conformal_requires_two_null_vectors = 2,
  ga_status_id_inverse_zero_denominator = 3,
  ga_status_id_accumulator_stack_exhausted = 4
};

typedef struct {
  ga_blade_id_t id;
  double coeff;
} ga_term_t;

typedef struct {
  ga_term_t *term;
  uint64_t term_count;
} ga_mv_t;

typedef struct {
  double coeff[ga_blade_count];
  ga_blade_id_t used_id[ga_blade_count];
  uint64_t used_count;
} ga_accumulator_t;

typedef struct {
  uint64_t dimensions;
  uint64_t is_conformal;
  uint64_t is_symmetric;
  uint64_t is_diagonal;
  uint64_t null_vector_count;
  uint64_t null_vector_start;
  ga_blade_id_t id_null;
  ga_blade_id_t pseudoscalar_id;
  ga_blade_id_t no_id;
  ga_blade_id_t ni_id;
  uint64_t no_bit_index;
  uint64_t ni_bit_index;
  uint64_t no_index;
  uint64_t ni_index;
  double metric[ga_max_dimensions][ga_max_dimensions];
  uint8_t grade_by_id[ga_blade_count];
  uint64_t accumulator_count;
  ga_accumulator_t accumulator[ga_accumulator_stack_size];
} ga_t;

char *ga_status_description(status_t status);

status_t ga_init_full(ga_t *ga, double metric[ga_max_dimensions][ga_max_dimensions], uint64_t dimensions, uint64_t is_conformal);
status_t ga_init_diagonal(ga_t *ga, double *diagonal_metric, uint64_t dimensions, uint64_t is_conformal);
void ga_uninit(ga_t *ga);

status_t ga_add(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out);
status_t ga_subtract(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out);

void ga_involute(ga_t *ga, ga_mv_t *in, ga_mv_t *out);
void ga_reverse(ga_t *ga, ga_mv_t *in, ga_mv_t *out);
void ga_conjugate(ga_t *ga, ga_mv_t *in, ga_mv_t *out);

status_t ga_ep_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out);
status_t ga_gp_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out);
status_t ga_ip_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out);
status_t ga_inverse(ga_t *ga, ga_mv_t *in, ga_mv_t *reverse_mv, ga_mv_t *denominator_mv, ga_mv_t *out);

void ga_mv_from_string(ga_t *ga, char *string, ga_mv_t *out);
void ga_mv_to_string(ga_t *ga, ga_mv_t *in, char *buffer);

void ga_no(ga_t *ga, double coeff, ga_mv_t *out);
void ga_ni(ga_t *ga, double coeff, ga_mv_t *out);
void ga_point(ga_t *ga, double *euclidean_coeff, ga_mv_t *out);
void ga_point_euclidean(ga_t *ga, ga_mv_t *in, double *euclidean_coeff);
void ga_rotor(ga_t *ga, double *coeff, ga_mv_t *out);

#endif