#include <sph-ga/sph-ga.h>
#include <stdio.h>
#include <stdlib.h>

char ga_status_group[] = "ga";

char *ga_status_description(status_t status)
{
  if (!status.id) {
    return "success";
  }
  if (status.group != ga_status_group) {
    return "";
  }
  if (status.id == ga_status_id_metric_not_symmetric) {
    return "metric not symmetric";
  }
  if (status.id == ga_status_id_conformal_requires_two_null_vectors) {
    return "conformal requires two null vectors";
  }
  if (status.id == ga_status_id_inverse_zero_denominator) {
    return "inverse zero denominator";
  }
  if (status.id == ga_status_id_accumulator_stack_exhausted) {
    return "accumulator stack exhausted";
  }
  return "";
}

double ga_determinant_generic(double matrix[ga_max_dimensions][ga_max_dimensions], uint64_t dimensions)
{
  double minor[ga_max_dimensions][ga_max_dimensions];
  double result;
  double sign;
  uint64_t column_index;
  uint64_t source_row;
  uint64_t source_column;
  uint64_t minor_column;

  if (dimensions == 1) {
    return matrix[0][0];
  }

  result = 0.0;
  column_index = 0;
  while (column_index < dimensions) {
    source_row = 1;
    while (source_row < dimensions) {
      source_column = 0;
      minor_column = 0;
      while (source_column < dimensions) {
        if (source_column != column_index) {
          minor[source_row - 1][minor_column] = matrix[source_row][source_column];
          minor_column += 1;
        }
        source_column += 1;
      }
      source_row += 1;
    }
    sign = (column_index & 1) ? -1.0 : 1.0;
    result += sign * matrix[0][column_index] * ga_determinant_generic(minor, dimensions - 1);
    column_index += 1;
  }

  return result;
}

double ga_determinant(double matrix[ga_max_dimensions][ga_max_dimensions], uint64_t dimensions)
{
  if (dimensions == 1) {
    return matrix[0][0];
  }
  if (dimensions == 2) {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  }
  if (dimensions == 3) {
    return
      matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
      matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
      matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
  }
  if (dimensions == 4) {
    return
      matrix[0][0] * (
        matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
        matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
        matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])
      ) -
      matrix[0][1] * (
        matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
        matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
        matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])
      ) +
      matrix[0][2] * (
        matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
        matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
        matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
      ) -
      matrix[0][3] * (
        matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
        matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
        matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
      );
  }
  return ga_determinant_generic(matrix, dimensions);
}

void ga_metric_properties(double metric[ga_max_dimensions][ga_max_dimensions], uint64_t dimensions, uint64_t *is_symmetric, uint64_t *is_diagonal, uint64_t *null_vector_count)
{
  uint64_t row_index;
  uint64_t column_index;
  uint64_t first_zero;
  double diagonal_value;

  *is_symmetric = 1;
  *is_diagonal = 1;
  *null_vector_count = 0;
  first_zero = dimensions;

  row_index = 0;
  while (row_index < dimensions) {
    diagonal_value = metric[row_index][row_index];
    if (!diagonal_value) {
      *null_vector_count += 1;
      if (first_zero == dimensions) {
        first_zero = row_index;
      }
    } else {
      if ((first_zero < dimensions) || !((diagonal_value == 1.0) || (diagonal_value == -1.0))) {
        *is_diagonal = 0;
        break;
      }
    }

    column_index = row_index + 1;
    while (column_index < dimensions) {
      if (row_index < first_zero) {
        if (metric[row_index][column_index] != metric[column_index][row_index]) {
          *is_symmetric = 0;
          *is_diagonal = 0;
          return;
        }
        if (metric[row_index][column_index] != 0.0) {
          *is_diagonal = 0;
          break;
        }
      }
      column_index += 1;
    }

    if (!*is_diagonal) {
      break;
    }
    row_index += 1;
  }

  row_index = first_zero;
  while (row_index < dimensions) {
    if (metric[row_index][row_index] != 0.0) {
      *is_symmetric = 0;
      *is_diagonal = 0;
      return;
    }
    row_index += 1;
  }
}

void ga_init_grade_table(ga_t *ga)
{
  ga_blade_id_t id;
  ga_blade_id_t value;
  uint64_t grade;
  uint64_t accumulator_index;

  id = 0;
  while (id < ga_blade_count) {
    value = id;
    grade = 0;
    while (value) {
      value &= value - 1;
      grade += 1;
    }
    ga->grade_by_id[id] = (uint8_t) grade;
    id += 1;
  }

  accumulator_index = 0;
  while (accumulator_index < ga_accumulator_stack_size) {
    id = 0;
    while (id < ga_blade_count) {
      ga->accumulator[accumulator_index].coeff[id] = 0.0;
      id += 1;
    }
    ga->accumulator[accumulator_index].used_count = 0;
    accumulator_index += 1;
  }
}

void ga_id_to_bit_indices(ga_blade_id_t id, uint64_t *bit_index, uint64_t *bit_index_count)
{
  uint64_t candidate;
  uint64_t count;

  candidate = 0;
  count = 0;
  while (candidate < ga_max_dimensions) {
    if (id & ((((ga_blade_id_t) 1)) << candidate)) {
      bit_index[count] = candidate;
      count += 1;
    }
    candidate += 1;
  }
  *bit_index_count = count;
}

uint64_t ga_id_to_indices(ga_blade_id_t id, uint64_t *index)
{
  uint64_t bit_index;
  uint64_t index_count;

  bit_index = 0;
  index_count = 0;
  while (bit_index < ga_max_dimensions) {
    if (id & ((((ga_blade_id_t) 1)) << bit_index)) {
      index[index_count] = bit_index + 1;
      index_count += 1;
    }
    bit_index += 1;
  }

  return index_count;
}

void ga_flat_metric_to_full(double *diagonal_metric, uint64_t dimensions, double full_metric[ga_max_dimensions][ga_max_dimensions])
{
  uint64_t row_index;
  uint64_t column_index;

  row_index = 0;
  while (row_index < dimensions) {
    column_index = 0;
    while (column_index < dimensions) {
      full_metric[row_index][column_index] = 0.0;
      column_index += 1;
    }
    full_metric[row_index][row_index] = diagonal_metric[row_index];
    row_index += 1;
  }
}

void ga_uninit(ga_t *ga)
{
  ga->dimensions = 0;
}

status_t ga_init_full(ga_t *ga, double metric[ga_max_dimensions][ga_max_dimensions], uint64_t dimensions, uint64_t is_conformal)
{
  status_declare;
  uint64_t row_index;
  uint64_t column_index;
  uint64_t base_dimensions;

  ga_uninit(ga);
  ga_init_grade_table(ga);

  row_index = 0;
  while (row_index < ga_max_dimensions) {
    column_index = 0;
    while (column_index < ga_max_dimensions) {
      ga->metric[row_index][column_index] = 0.0;
      column_index += 1;
    }
    row_index += 1;
  }

  ga->is_conformal = 0;
  ga->is_symmetric = 0;
  ga->is_diagonal = 0;
  ga->null_vector_count = 0;
  ga->null_vector_start = 0;
  ga->id_null = 0;
  ga->pseudoscalar_id = 0;
  ga->no_id = 0;
  ga->ni_id = 0;
  ga->no_bit_index = 0;
  ga->ni_bit_index = 0;
  ga->no_index = 0;
  ga->ni_index = 0;
  ga->accumulator_count = 0;

  base_dimensions = dimensions;
  if (is_conformal) {
    dimensions += 2;
  }

  row_index = 0;
  while (row_index < base_dimensions) {
    column_index = 0;
    while (column_index < base_dimensions) {
      ga->metric[row_index][column_index] = metric[row_index][column_index];
      column_index += 1;
    }
    row_index += 1;
  }

  if (is_conformal) {
    ga->metric[dimensions - 2][dimensions - 1] = -1.0;
    ga->metric[dimensions - 1][dimensions - 2] = -1.0;
  }

  ga_metric_properties(ga->metric, dimensions, &ga->is_symmetric, &ga->is_diagonal, &ga->null_vector_count);

  if (!ga->is_symmetric) {
    status_set_goto(ga_status_group, ga_status_id_metric_not_symmetric);
  }
  if (is_conformal) {
    if (ga->null_vector_count != 2) {
      status_set_goto(ga_status_group, ga_status_id_conformal_requires_two_null_vectors);
    }
  }

  ga->dimensions = dimensions;
  ga->is_conformal = is_conformal;
  ga->pseudoscalar_id = ((((ga_blade_id_t) 1)) << dimensions) - 1;

  if (ga->null_vector_count) {
    ga->null_vector_start = dimensions - ga->null_vector_count;
    ga->id_null = (((((ga_blade_id_t) 1)) << ga->null_vector_count) - 1) << ga->null_vector_start;
  }

  if (is_conformal) {
    ga->no_bit_index = dimensions - 2;
    ga->ni_bit_index = dimensions - 1;
    ga->no_id = (((ga_blade_id_t) 1)) << ga->no_bit_index;
    ga->ni_id = (((ga_blade_id_t) 1)) << ga->ni_bit_index;
    ga->no_index = ga->no_bit_index + 1;
    ga->ni_index = ga->ni_bit_index + 1;
  }

exit:
  if (status_is_failure) {
    ga_uninit(ga);
  }
  status_return;
}

status_t ga_init_diagonal(ga_t *ga, double *diagonal_metric, uint64_t dimensions, uint64_t is_conformal)
{
  status_declare;
  double full_metric[ga_max_dimensions][ga_max_dimensions];

  ga_flat_metric_to_full(diagonal_metric, dimensions, full_metric);
  status_require((ga_init_full(ga, full_metric, dimensions, is_conformal)));

exit:
  status_return;
}

void ga_involute(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  uint64_t term_index;
  uint64_t grade;
  double factor;

  term_index = 0;
  while (term_index < in->term_count) {
    grade = ga->grade_by_id[in->term[term_index].id];
    factor = (grade & 1) ? -1.0 : 1.0;
    out->term[term_index].id = in->term[term_index].id;
    out->term[term_index].coeff = in->term[term_index].coeff * factor;
    term_index += 1;
  }
  out->term_count = in->term_count;
}

void ga_reverse(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  uint64_t term_index;
  uint64_t grade;
  double factor;

  term_index = 0;
  while (term_index < in->term_count) {
    grade = ga->grade_by_id[in->term[term_index].id];
    factor = ((((grade * (grade - 1)) >> 1) & 1) ? -1.0 : 1.0);
    out->term[term_index].id = in->term[term_index].id;
    out->term[term_index].coeff = in->term[term_index].coeff * factor;
    term_index += 1;
  }
  out->term_count = in->term_count;
}

void ga_conjugate(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  uint64_t term_index;
  uint64_t grade;
  double factor;

  term_index = 0;
  while (term_index < in->term_count) {
    grade = ga->grade_by_id[in->term[term_index].id];
    factor = ((((grade * (grade + 1)) >> 1) & 1) ? -1.0 : 1.0);
    out->term[term_index].id = in->term[term_index].id;
    out->term[term_index].coeff = in->term[term_index].coeff * factor;
    term_index += 1;
  }
  out->term_count = in->term_count;
}

void ga_array_remove_pair(uint64_t *value, uint64_t *value_count, uint64_t left_index, uint64_t right_index)
{
  uint64_t index;

  if (left_index < right_index) {
    index = right_index;
    while (index + 1 < *value_count) {
      value[index] = value[index + 1];
      index += 1;
    }
    *value_count -= 1;

    index = left_index;
    while (index + 1 < *value_count) {
      value[index] = value[index + 1];
      index += 1;
    }
    *value_count -= 1;
  } else {
    index = left_index;
    while (index + 1 < *value_count) {
      value[index] = value[index + 1];
      index += 1;
    }
    *value_count -= 1;

    index = right_index;
    while (index + 1 < *value_count) {
      value[index] = value[index + 1];
      index += 1;
    }
    *value_count -= 1;
  }
}

double ga_sign_indices(uint64_t *index, uint64_t index_count)
{
  uint64_t left_index;
  uint64_t right_index;
  uint64_t inversion_count;

  inversion_count = 0;
  left_index = 0;
  while (left_index < index_count) {
    right_index = left_index + 1;
    while (right_index < index_count) {
      if (index[left_index] > index[right_index]) {
        inversion_count += 1;
      }
      right_index += 1;
    }
    left_index += 1;
  }

  return (inversion_count & 1) ? -1.0 : 1.0;
}

double ga_sign_sorted(uint64_t *left, uint64_t left_count, uint64_t *right, uint64_t right_count)
{
  uint64_t left_index;
  uint64_t right_index;
  uint64_t inversion_count;

  left_index = 0;
  right_index = 0;
  inversion_count = 0;

  while ((left_index < left_count) && (right_index < right_count)) {
    if (left[left_index] <= right[right_index]) {
      left_index += 1;
    } else {
      inversion_count += 1;
      right_index += 1;
    }
  }

  return (inversion_count & 1) ? -1.0 : 1.0;
}

status_t ga_add(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out)
{
  status_declare;
  ga_accumulator_t *accumulator;
  uint64_t term_index;
  uint64_t used_index;
  ga_blade_id_t id;
  double coeff;

  if (ga->accumulator_count >= ga_accumulator_stack_size) {
    status_set_goto(ga_status_group, ga_status_id_accumulator_stack_exhausted);
  }

  accumulator = ga->accumulator + ga->accumulator_count;
  ga->accumulator_count += 1;
  used_index = 0;
  while (used_index < accumulator->used_count) {
    accumulator->coeff[accumulator->used_id[used_index]] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;

  term_index = 0;
  while (term_index < left->term_count) {
    id = left->term[term_index].id;
    coeff = left->term[term_index].coeff;
    if (!accumulator->coeff[id]) {
      accumulator->used_id[accumulator->used_count] = id;
      accumulator->used_count += 1;
    }
    accumulator->coeff[id] += coeff;
    term_index += 1;
  }

  term_index = 0;
  while (term_index < right->term_count) {
    id = right->term[term_index].id;
    coeff = right->term[term_index].coeff;
    if (!accumulator->coeff[id]) {
      accumulator->used_id[accumulator->used_count] = id;
      accumulator->used_count += 1;
    }
    accumulator->coeff[id] += coeff;
    term_index += 1;
  }

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      out->term[term_index].id = id;
      out->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  out->term_count = term_index;

  ga->accumulator_count -= 1;

exit:
  status_return;
}

status_t ga_subtract(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out)
{
  status_declare;
  ga_accumulator_t *accumulator;
  uint64_t term_index;
  uint64_t used_index;
  ga_blade_id_t id;
  double coeff;

  if (ga->accumulator_count >= ga_accumulator_stack_size) {
    status_set_goto(ga_status_group, ga_status_id_accumulator_stack_exhausted);
  }

  accumulator = ga->accumulator + ga->accumulator_count;
  ga->accumulator_count += 1;
  used_index = 0;
  while (used_index < accumulator->used_count) {
    accumulator->coeff[accumulator->used_id[used_index]] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;

  term_index = 0;
  while (term_index < left->term_count) {
    id = left->term[term_index].id;
    coeff = left->term[term_index].coeff;
    if (!accumulator->coeff[id]) {
      accumulator->used_id[accumulator->used_count] = id;
      accumulator->used_count += 1;
    }
    accumulator->coeff[id] += coeff;
    term_index += 1;
  }

  term_index = 0;
  while (term_index < right->term_count) {
    id = right->term[term_index].id;
    coeff = -right->term[term_index].coeff;
    if (!accumulator->coeff[id]) {
      accumulator->used_id[accumulator->used_count] = id;
      accumulator->used_count += 1;
    }
    accumulator->coeff[id] += coeff;
    term_index += 1;
  }

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      out->term[term_index].id = id;
      out->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  out->term_count = term_index;

  ga->accumulator_count -= 1;

exit:
  status_return;
}

status_t ga_ep_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out)
{
  status_declare;
  ga_accumulator_t *accumulator;
  uint64_t left_term_index;
  uint64_t right_term_index;
  ga_blade_id_t left_id;
  ga_blade_id_t right_id;
  ga_blade_id_t id;
  uint64_t left_bit_index[ga_max_dimensions];
  uint64_t right_bit_index[ga_max_dimensions];
  uint64_t left_bit_index_count;
  uint64_t right_bit_index_count;
  uint64_t used_index;
  uint64_t term_index;
  double sign;
  double coeff;

  if (ga->accumulator_count >= ga_accumulator_stack_size) {
    status_set_goto(ga_status_group, ga_status_id_accumulator_stack_exhausted);
  }

  accumulator = ga->accumulator + ga->accumulator_count;
  ga->accumulator_count += 1;
  used_index = 0;
  while (used_index < accumulator->used_count) {
    accumulator->coeff[accumulator->used_id[used_index]] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    left_id = left->term[left_term_index].id;
    left_bit_index_count = 0;
    if (left_id) {
      ga_id_to_bit_indices(left_id, left_bit_index, &left_bit_index_count);
    }

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      right_id = right->term[right_term_index].id;
      if (!(left_id & right_id)) {
        id = left_id | right_id;
        if (id) {
          sign = 1.0;
          if (right_id) {
            ga_id_to_bit_indices(right_id, right_bit_index, &right_bit_index_count);
            sign = ga_sign_sorted(left_bit_index, left_bit_index_count, right_bit_index, right_bit_index_count);
          }
          coeff = sign * left->term[left_term_index].coeff * right->term[right_term_index].coeff;
          if (coeff) {
            if (!accumulator->coeff[id]) {
              accumulator->used_id[accumulator->used_count] = id;
              accumulator->used_count += 1;
            }
            accumulator->coeff[id] += coeff;
          }
        }
      }
      right_term_index += 1;
    }

    left_term_index += 1;
  }

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      out->term[term_index].id = id;
      out->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  out->term_count = term_index;

  ga->accumulator_count -= 1;

exit:
  status_return;
}

status_t ga_gp_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out)
{
  status_declare;
  ga_accumulator_t *accumulator;
  uint64_t left_term_index;
  uint64_t right_term_index;
  uint64_t left_bit_index[ga_max_dimensions];
  uint64_t right_bit_index[ga_max_dimensions];
  uint64_t merged_bit_index[ga_max_dimensions * 2];
  uint64_t left_bit_index_count;
  uint64_t right_bit_index_count;
  uint64_t merged_bit_index_count;
  uint64_t index;
  uint64_t row_index;
  uint64_t column_index;
  uint64_t used_index;
  uint64_t term_index;
  uint64_t changed;
  ga_blade_id_t id;
  double sign;
  double factor;
  double metric_value;
  double coeff;

  if (ga->accumulator_count >= ga_accumulator_stack_size) {
    status_set_goto(ga_status_group, ga_status_id_accumulator_stack_exhausted);
  }

  accumulator = ga->accumulator + ga->accumulator_count;
  ga->accumulator_count += 1;
  used_index = 0;
  while (used_index < accumulator->used_count) {
    accumulator->coeff[accumulator->used_id[used_index]] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    ga_id_to_bit_indices(left->term[left_term_index].id, left_bit_index, &left_bit_index_count);

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      ga_id_to_bit_indices(right->term[right_term_index].id, right_bit_index, &right_bit_index_count);

      index = 0;
      while (index < left_bit_index_count) {
        merged_bit_index[index] = left_bit_index[index];
        index += 1;
      }
      merged_bit_index_count = left_bit_index_count;
      index = 0;
      while (index < right_bit_index_count) {
        merged_bit_index[merged_bit_index_count] = right_bit_index[index];
        merged_bit_index_count += 1;
        index += 1;
      }

      sign = ga_sign_indices(merged_bit_index, merged_bit_index_count);
      factor = 1.0;

      changed = 1;
      while (changed) {
        changed = 0;
        row_index = 0;
        while (row_index < merged_bit_index_count) {
          column_index = row_index + 1;
          while (column_index < merged_bit_index_count) {
            if (merged_bit_index[row_index] != merged_bit_index[column_index]) {
              metric_value = ga->metric[merged_bit_index[row_index]][merged_bit_index[column_index]];
              if (metric_value != 0.0) {
                factor *= metric_value;
                ga_array_remove_pair(merged_bit_index, &merged_bit_index_count, row_index, column_index);
                changed = 1;
                break;
              }
            }
            column_index += 1;
          }
          if (changed) {
            break;
          }
          row_index += 1;
        }
      }

      changed = 1;
      while (changed) {
        changed = 0;
        row_index = 0;
        while (row_index < merged_bit_index_count) {
          column_index = row_index + 1;
          while (column_index < merged_bit_index_count) {
            if (merged_bit_index[row_index] == merged_bit_index[column_index]) {
              metric_value = ga->metric[merged_bit_index[row_index]][merged_bit_index[row_index]];
              factor *= metric_value;
              ga_array_remove_pair(merged_bit_index, &merged_bit_index_count, row_index, column_index);
              changed = 1;
              break;
            }
            column_index += 1;
          }
          if (changed) {
            break;
          }
          row_index += 1;
        }
      }

      coeff = left->term[left_term_index].coeff * right->term[right_term_index].coeff * sign * factor;
      if (coeff) {
        id = 0;
        index = 0;
        while (index < merged_bit_index_count) {
          id |= (((ga_blade_id_t) 1)) << merged_bit_index[index];
          index += 1;
        }
        if (!accumulator->coeff[id]) {
          accumulator->used_id[accumulator->used_count] = id;
          accumulator->used_count += 1;
        }
        accumulator->coeff[id] += coeff;
      }

      right_term_index += 1;
    }

    left_term_index += 1;
  }

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      out->term[term_index].id = id;
      out->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  out->term_count = term_index;

  ga->accumulator_count -= 1;

exit:
  status_return;
}

double ga_sign_of_choice(uint64_t *choice, uint64_t choice_count)
{
  uint64_t index;
  uint64_t sum;

  sum = 0;
  index = 0;
  while (index < choice_count) {
    sum += choice[index] - index;
    index += 1;
  }
  return (sum & 1) ? -1.0 : 1.0;
}

void ga_choose_rec(uint64_t *space, uint64_t space_count, uint64_t choose_count, uint64_t space_index, uint64_t chosen_count, uint64_t *chosen, uint64_t *combo_count, uint64_t combo[ga_blade_count][ga_max_dimensions])
{
  uint64_t copy_index;

  if (chosen_count == choose_count) {
    copy_index = 0;
    while (copy_index < choose_count) {
      combo[*combo_count][copy_index] = chosen[copy_index];
      copy_index += 1;
    }
    *combo_count += 1;
    return;
  }

  if (space_count - space_index < choose_count - chosen_count) {
    return;
  }

  chosen[chosen_count] = space[space_index];
  ga_choose_rec(space, space_count, choose_count, space_index + 1, chosen_count + 1, chosen, combo_count, combo);
  ga_choose_rec(space, space_count, choose_count, space_index + 1, chosen_count, chosen, combo_count, combo);
}

status_t ga_ip_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_mv_t *out)
{
  status_declare;
  ga_accumulator_t *accumulator;
  uint64_t left_term_index;
  uint64_t right_term_index;
  uint64_t left_bit_index[ga_max_dimensions];
  uint64_t right_bit_index[ga_max_dimensions];
  uint64_t left_bit_index_count;
  uint64_t right_bit_index_count;
  uint64_t space[ga_max_dimensions];
  uint64_t chosen[ga_max_dimensions];
  uint64_t combo[ga_blade_count][ga_max_dimensions];
  uint64_t combo_count;
  uint64_t combo_index;
  uint64_t combo_position;
  uint64_t remainder_bit_index[ga_max_dimensions];
  uint64_t remainder_bit_index_count;
  uint64_t index;
  uint64_t found;
  uint64_t row_index;
  uint64_t column_index;
  uint64_t used_index;
  uint64_t term_index;
  uint64_t sign_power;
  ga_blade_id_t id;
  double sign;
  double coeff;
  double matrix[ga_max_dimensions][ga_max_dimensions];
  double scalar;

  if (ga->accumulator_count >= ga_accumulator_stack_size) {
    status_set_goto(ga_status_group, ga_status_id_accumulator_stack_exhausted);
  }

  accumulator = ga->accumulator + ga->accumulator_count;
  ga->accumulator_count += 1;
  used_index = 0;
  while (used_index < accumulator->used_count) {
    accumulator->coeff[accumulator->used_id[used_index]] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    ga_id_to_bit_indices(left->term[left_term_index].id, left_bit_index, &left_bit_index_count);

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      ga_id_to_bit_indices(right->term[right_term_index].id, right_bit_index, &right_bit_index_count);

      if (right_bit_index_count >= left_bit_index_count) {
        if (right_bit_index_count == left_bit_index_count) {
          if (left_bit_index_count == 1) {
            coeff =
              left->term[left_term_index].coeff *
              right->term[right_term_index].coeff *
              ga->metric[left_bit_index[0]][right_bit_index[0]];
            if (coeff) {
              if (!accumulator->coeff[0]) {
                accumulator->used_id[accumulator->used_count] = 0;
                accumulator->used_count += 1;
              }
              accumulator->coeff[0] += coeff;
            }
          } else {
            row_index = 0;
            while (row_index < left_bit_index_count) {
              column_index = 0;
              while (column_index < left_bit_index_count) {
                matrix[row_index][column_index] = ga->metric[left_bit_index[row_index]][right_bit_index[column_index]];
                column_index += 1;
              }
              row_index += 1;
            }

            scalar = ga_determinant(matrix, left_bit_index_count);
            sign_power = (left_bit_index_count * (left_bit_index_count - 1)) / 2;
            if (sign_power & 1) {
              scalar = -scalar;
            }
            coeff = -scalar * left->term[left_term_index].coeff * right->term[right_term_index].coeff;

            if (coeff) {
              if (!accumulator->coeff[0]) {
                accumulator->used_id[accumulator->used_count] = 0;
                accumulator->used_count += 1;
              }
              accumulator->coeff[0] += coeff;
            }
          }
        } else {
          index = 0;
          while (index < right_bit_index_count) {
            space[index] = index;
            index += 1;
          }

          combo_count = 0;
          ga_choose_rec(space, right_bit_index_count, left_bit_index_count, 0, 0, chosen, &combo_count, combo);

          combo_index = 0;
          while (combo_index < combo_count) {
            sign = ga_sign_of_choice(combo[combo_index], left_bit_index_count);
            coeff = 1.0;

            index = 0;
            while (index < left_bit_index_count) {
              coeff *= ga->metric[left_bit_index[index]][right_bit_index[combo[combo_index][index]]];
              index += 1;
            }

            remainder_bit_index_count = 0;
            index = 0;
            while (index < right_bit_index_count) {
              combo_position = 0;
              found = 0;
              while (combo_position < left_bit_index_count) {
                if (combo[combo_index][combo_position] == index) {
                  found = 1;
                  break;
                }
                combo_position += 1;
              }
              if (!found) {
                remainder_bit_index[remainder_bit_index_count] = right_bit_index[index];
                remainder_bit_index_count += 1;
              }
              index += 1;
            }

            id = 0;
            index = 0;
            while (index < remainder_bit_index_count) {
              id |= (((ga_blade_id_t) 1)) << remainder_bit_index[index];
              index += 1;
            }

            coeff *= left->term[left_term_index].coeff * right->term[right_term_index].coeff * sign;
            if (coeff) {
              if (!accumulator->coeff[id]) {
                accumulator->used_id[accumulator->used_count] = id;
                accumulator->used_count += 1;
              }
              accumulator->coeff[id] += coeff;
            }

            combo_index += 1;
          }
        }
      }

      right_term_index += 1;
    }

    left_term_index += 1;
  }

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      out->term[term_index].id = id;
      out->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  out->term_count = term_index;

  ga->accumulator_count -= 1;

exit:
  status_return;
}

status_t ga_inverse(ga_t *ga, ga_mv_t *in, ga_mv_t *reverse_mv, ga_mv_t *denominator_mv, ga_mv_t *out)
{
  status_declare;
  uint64_t term_index;
  double denominator;

  denominator = 0.0;

  ga_reverse(ga, in, reverse_mv);
  status_require((ga_gp_one(ga, in, reverse_mv, denominator_mv)));

  term_index = 0;
  while (term_index < denominator_mv->term_count) {
    if (!denominator_mv->term[term_index].id) {
      denominator = denominator_mv->term[term_index].coeff;
      break;
    }
    term_index += 1;
  }

  if (!denominator) {
    status_set_goto(ga_status_group, ga_status_id_inverse_zero_denominator);
  }

  term_index = 0;
  while (term_index < reverse_mv->term_count) {
    out->term[term_index].id = reverse_mv->term[term_index].id;
    out->term[term_index].coeff = reverse_mv->term[term_index].coeff / denominator;
    term_index += 1;
  }
  out->term_count = reverse_mv->term_count;

exit:
  status_return;
}

#define ga_skip_space(string) while (**(string) == ' ') { *(string) += 1; }

void ga_mv_from_string(ga_t *ga, char *string, ga_mv_t *out)
{
  char *read;
  uint64_t term_index;
  uint64_t out_index;
  uint64_t index[ga_max_dimensions];
  uint64_t index_count;
  uint64_t value;
  uint64_t digit_count;
  ga_blade_id_t id;
  double coeff;
  double scale;
  double whole_part;
  double fractional_part;
  int sign;
  int found;

  read = string;
  out->term_count = 0;

  while (1) {
    ga_skip_space((&read));
    if (!*read) {
      break;
    }

    sign = 1;
    if (*read == '+') {
      read += 1;
      ga_skip_space((&read));
    }
    if (*read == '-') {
      sign = -1;
      read += 1;
      ga_skip_space((&read));
    }

    coeff = 1.0;
    whole_part = 0.0;
    fractional_part = 0.0;
    scale = 1.0;
    digit_count = 0;

    while ((*read >= '0') && (*read <= '9')) {
      whole_part *= 10.0;
      whole_part += *read - '0';
      read += 1;
      digit_count += 1;
    }

    if (*read == '.') {
      read += 1;
      while ((*read >= '0') && (*read <= '9')) {
        fractional_part *= 10.0;
        fractional_part += *read - '0';
        scale *= 10.0;
        read += 1;
        digit_count += 1;
      }
    }

    if (digit_count) {
      coeff = whole_part + (fractional_part / scale);
    }

    id = 0;
    if (*read == 'e') {
      read += 1;
      index_count = 0;
      while (1) {
        value = 0;
        while ((*read >= '0') && (*read <= '9')) {
          value *= 10;
          value += *read - '0';
          read += 1;
        }
        index[index_count] = value;
        index_count += 1;
        if (*read != '_') {
          break;
        }
        read += 1;
      }

      term_index = 0;
      while (term_index < index_count) {
        id |= (((ga_blade_id_t) 1)) << (index[term_index] - 1);
        term_index += 1;
      }
    } else if ((read[0] == 'n') && (read[1] == 'o')) {
      id = ga->no_id;
      read += 2;
    } else if ((read[0] == 'n') && (read[1] == 'i')) {
      id = ga->ni_id;
      read += 2;
    }

    coeff *= sign;

    found = 0;
    out_index = 0;
    while (out_index < out->term_count) {
      if (out->term[out_index].id == id) {
        out->term[out_index].coeff += coeff;
        found = 1;
        break;
      }
      out_index += 1;
    }

    if (!found && coeff) {
      out->term[out->term_count].id = id;
      out->term[out->term_count].coeff = coeff;
      out->term_count += 1;
    }

    ga_skip_space((&read));
    if (!*read) {
      break;
    }
  }

  out_index = 0;
  term_index = 0;
  while (term_index < out->term_count) {
    if (out->term[term_index].coeff) {
      if (out_index != term_index) {
        out->term[out_index] = out->term[term_index];
      }
      out_index += 1;
    }
    term_index += 1;
  }
  out->term_count = out_index;
}

void ga_mv_to_string(ga_t *ga, ga_mv_t *in, char *buffer)
{
  char *write;
  uint64_t first;
  uint64_t search_index;
  uint64_t index[ga_max_dimensions];
  uint64_t index_count;
  uint64_t index_position;
  ga_blade_id_t id;
  ga_blade_id_t limit;
  double coeff;

  write = buffer;
  first = 1;
  limit = (((ga_blade_id_t) 1)) << ga->dimensions;

  id = 0;
  while (id < limit) {
    coeff = 0.0;
    search_index = 0;
    while (search_index < in->term_count) {
      if (in->term[search_index].id == id) {
        coeff = in->term[search_index].coeff;
        break;
      }
      search_index += 1;
    }

    if (coeff) {
      if (!first) {
        write += sprintf(write, " + ");
      }

      if (!id) {
        write += sprintf(write, "%g", coeff);
      } else {
        if (coeff == -1.0) {
          write += sprintf(write, "-");
        } else if (coeff != 1.0) {
          write += sprintf(write, "%g", coeff);
        }

        index_count = ga_id_to_indices(id, index);
        write += sprintf(write, "e%" PRIu64, index[0]);
        index_position = 1;
        while (index_position < index_count) {
          write += sprintf(write, "_%" PRIu64, index[index_position]);
          index_position += 1;
        }
      }

      first = 0;
    }

    id += 1;
  }

  if (first) {
    sprintf(write, "0");
  } else {
    *write = 0;
  }
}

void ga_no(ga_t *ga, double coeff, ga_mv_t *out)
{
  if (!coeff) {
    out->term_count = 0;
    return;
  }

  out->term[0].id = ga->no_id;
  out->term[0].coeff = coeff;
  out->term_count = 1;
}

void ga_ni(ga_t *ga, double coeff, ga_mv_t *out)
{
  if (!coeff) {
    out->term_count = 0;
    return;
  }

  out->term[0].id = ga->ni_id;
  out->term[0].coeff = coeff;
  out->term_count = 1;
}

void ga_point(ga_t *ga, double *euclidean_coeff, ga_mv_t *out)
{
  uint64_t euclidean_dimensions;
  uint64_t dimension_index;
  uint64_t term_index;
  double ni_coeff;

  euclidean_dimensions = ga->dimensions - 2;
  ni_coeff = 0.0;
  term_index = 0;

  dimension_index = 0;
  while (dimension_index < euclidean_dimensions) {
    if (euclidean_coeff[dimension_index]) {
      out->term[term_index].id = (((ga_blade_id_t) 1)) << dimension_index;
      out->term[term_index].coeff = euclidean_coeff[dimension_index];
      term_index += 1;
      ni_coeff += euclidean_coeff[dimension_index] * euclidean_coeff[dimension_index];
    }
    dimension_index += 1;
  }

  out->term[term_index].id = ga->no_id;
  out->term[term_index].coeff = 1.0;
  term_index += 1;

  ni_coeff *= 0.5;
  if (ni_coeff) {
    out->term[term_index].id = ga->ni_id;
    out->term[term_index].coeff = ni_coeff;
    term_index += 1;
  }

  out->term_count = term_index;
}

void ga_point_euclidean(ga_t *ga, ga_mv_t *in, double *euclidean_coeff)
{
  uint64_t euclidean_dimensions;
  uint64_t dimension_index;
  uint64_t term_index;
  double no_coeff;
  ga_blade_id_t id;

  euclidean_dimensions = ga->dimensions - 2;
  no_coeff = 0.0;

  dimension_index = 0;
  while (dimension_index < euclidean_dimensions) {
    euclidean_coeff[dimension_index] = 0.0;
    dimension_index += 1;
  }

  term_index = 0;
  while (term_index < in->term_count) {
    if (in->term[term_index].id == ga->no_id) {
      no_coeff = in->term[term_index].coeff;
      break;
    }
    term_index += 1;
  }

  term_index = 0;
  while (term_index < in->term_count) {
    id = in->term[term_index].id;
    dimension_index = 0;
    while (dimension_index < euclidean_dimensions) {
      if (id == ((((ga_blade_id_t) 1)) << dimension_index)) {
        euclidean_coeff[dimension_index] = in->term[term_index].coeff / no_coeff;
        break;
      }
      dimension_index += 1;
    }
    term_index += 1;
  }
}

void ga_rotor(ga_t *ga, double *coeff, ga_mv_t *out)
{
  uint64_t euclidean_dimensions;
  uint64_t left_index;
  uint64_t right_index;
  uint64_t coeff_index;
  uint64_t term_index;
  ga_blade_id_t id;

  euclidean_dimensions = ga->dimensions - 2;
  coeff_index = 0;
  term_index = 0;

  if (coeff[coeff_index]) {
    out->term[term_index].id = 0;
    out->term[term_index].coeff = coeff[coeff_index];
    term_index += 1;
  }
  coeff_index += 1;

  left_index = 0;
  while (left_index < euclidean_dimensions) {
    right_index = left_index + 1;
    while (right_index < euclidean_dimensions) {
      if (coeff[coeff_index]) {
        id = 0;
        id |= (((ga_blade_id_t) 1)) << left_index;
        id |= (((ga_blade_id_t) 1)) << right_index;
        out->term[term_index].id = id;
        out->term[term_index].coeff = coeff[coeff_index];
        term_index += 1;
      }
      coeff_index += 1;
      right_index += 1;
    }
    left_index += 1;
  }

  out->term_count = term_index;
}