#include <sph-ga/sph-ga.h>

char ga_status_group[] = "sph-ga";

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
  if (status.id == ga_status_id_conformal_requires_diagonal_metric) {
    return "conformal requires diagonal metric";
  }
  if (status.id == ga_status_id_inverse_zero_denominator) {
    return "inverse zero denominator";
  }
  return "";
}

void ga_uninit(ga_t *ga)
{
  ga->dimensions = 0;
}

void ga_accumulator_clear(ga_accumulator_t *accumulator)
{
  int64_t used_index;
  ga_blade_id_t id;

  used_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    accumulator->coeff[id] = 0.0;
    used_index += 1;
  }
  accumulator->used_count = 0;
}

void ga_accumulator_add(ga_accumulator_t *accumulator, ga_blade_id_t id, double coeff)
{
  if (!coeff) {
    return;
  }
  if (!accumulator->coeff[id]) {
    accumulator->used_id[accumulator->used_count] = id;
    accumulator->used_count += 1;
  }
  accumulator->coeff[id] += coeff;
}

void ga_accumulator_to_mv(ga_t *ga, ga_accumulator_t *accumulator, ga_mv_t *mv)
{
  int64_t used_index;
  int64_t term_index;
  ga_blade_id_t id;
  double coeff;

  (void) ga;

  used_index = 0;
  term_index = 0;
  while (used_index < accumulator->used_count) {
    id = accumulator->used_id[used_index];
    coeff = accumulator->coeff[id];
    if (coeff) {
      mv->term[term_index].id = id;
      mv->term[term_index].coeff = coeff;
      term_index += 1;
    }
    used_index += 1;
  }
  mv->term_count = term_index;
}

void ga_zero(ga_mv_t *mv)
{
  mv->term_count = 0;
}

void ga_scalar(double coeff, ga_mv_t *mv)
{
  if (coeff) {
    mv->term[0].id = 0;
    mv->term[0].coeff = coeff;
    mv->term_count = 1;
  } else {
    mv->term_count = 0;
  }
}

void ga_basis_blade(int64_t index, double coeff, ga_mv_t *mv)
{
  if (!coeff) {
    mv->term_count = 0;
    return;
  }
  if (!index) {
    mv->term[0].id = 0;
    mv->term[0].coeff = coeff;
    mv->term_count = 1;
    return;
  }
  mv->term[0].id = (((ga_blade_id_t) 1)) << (index - 1);
  mv->term[0].coeff = coeff;
  mv->term_count = 1;
}

void ga_scale(ga_mv_t *in, double scalar, ga_mv_t *out)
{
  int64_t term_index;

  term_index = 0;
  while (term_index < in->term_count) {
    out->term[term_index].id = in->term[term_index].id;
    out->term[term_index].coeff = in->term[term_index].coeff * scalar;
    term_index += 1;
  }
  out->term_count = in->term_count;
}

void ga_init_grade_table(ga_t *ga)
{
  ga_blade_id_t id;
  ga_blade_id_t value;
  int64_t grade;

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
}

ga_blade_id_t ga_id_from_bit_indices(int64_t *bit_index, int64_t bit_index_count)
{
  ga_blade_id_t id;
  int64_t bit_position;

  id = 0;
  bit_position = 0;
  while (bit_position < bit_index_count) {
    id |= (((ga_blade_id_t) 1)) << bit_index[bit_position];
    bit_position += 1;
  }
  return id;
}

ga_blade_id_t ga_id_from_indices(int64_t *index, int64_t index_count)
{
  ga_blade_id_t id;
  int64_t position;

  id = 0;
  position = 0;
  while (position < index_count) {
    if (index[position]) {
      id |= (((ga_blade_id_t) 1)) << (index[position] - 1);
    }
    position += 1;
  }
  return id;
}

void ga_id_to_bit_indices(ga_t *ga, ga_blade_id_t id, int64_t *bit_index, int64_t *bit_index_count)
{
  int64_t candidate;
  int64_t count;

  (void) ga;

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

void ga_flat_metric_to_full(double *diagonal_metric, int64_t dimensions, double full_metric[ga_max_dimensions][ga_max_dimensions])
{
  int64_t row_index;
  int64_t column_index;

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

void ga_metric_properties(double metric[ga_max_dimensions][ga_max_dimensions], int64_t dimensions, int64_t *is_symmetric, int64_t *is_diagonal, int64_t *null_vector_count)
{
  int64_t row_index;
  int64_t column_index;
  int64_t first_zero;
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

double ga_determinant_generic(double matrix[ga_max_dimensions][ga_max_dimensions], int64_t dimensions)
{
  double minor[ga_max_dimensions][ga_max_dimensions];
  double result;
  double sign;
  int64_t column_index;
  int64_t source_row;
  int64_t source_column;
  int64_t minor_column;

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

double ga_determinant(double matrix[ga_max_dimensions][ga_max_dimensions], int64_t dimensions)
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

void ga_array_remove_pair(int64_t *value, int64_t *value_count, int64_t left_index, int64_t right_index)
{
  int64_t index;

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

double ga_sign_indices(int64_t *index, int64_t index_count)
{
  int64_t left_index;
  int64_t right_index;
  int64_t inversion_count;

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

double ga_sign_sorted(int64_t *left, int64_t left_count, int64_t *right, int64_t right_count)
{
  int64_t left_index;
  int64_t right_index;
  int64_t inversion_count;

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

void ga_combine(ga_t *ga, ga_mv_t *left, ga_mv_t *right, double right_scale, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  int64_t term_index;

  ga_accumulator_clear(accumulator);

  term_index = 0;
  while (term_index < left->term_count) {
    ga_accumulator_add(accumulator, left->term[term_index].id, left->term[term_index].coeff);
    term_index += 1;
  }

  term_index = 0;
  while (term_index < right->term_count) {
    ga_accumulator_add(accumulator, right->term[term_index].id, right->term[term_index].coeff * right_scale);
    term_index += 1;
  }

  ga_accumulator_to_mv(ga, accumulator, out);
}

void ga_add(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  ga_combine(ga, left, right, 1.0, accumulator, out);
}

void ga_subtract(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  ga_combine(ga, left, right, -1.0, accumulator, out);
}

void ga_map_grade_factor(ga_t *ga, ga_mv_t *in, ga_mv_t *out, int64_t mode)
{
  int64_t term_index;
  int64_t grade;
  double factor;

  term_index = 0;
  while (term_index < in->term_count) {
    grade = ga->grade_by_id[in->term[term_index].id];
    if (!mode) {
      factor = (grade & 1) ? -1.0 : 1.0;
    } else if (mode == 1) {
      factor = ((((grade * (grade - 1)) >> 1) & 1) ? -1.0 : 1.0);
    } else {
      factor = ((((grade * (grade + 1)) >> 1) & 1) ? -1.0 : 1.0);
    }
    out->term[term_index].id = in->term[term_index].id;
    out->term[term_index].coeff = in->term[term_index].coeff * factor;
    term_index += 1;
  }
  out->term_count = in->term_count;
}

void ga_involute(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  ga_map_grade_factor(ga, in, out, 0);
}

void ga_reverse(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  ga_map_grade_factor(ga, in, out, 1);
}

void ga_conjugate(ga_t *ga, ga_mv_t *in, ga_mv_t *out)
{
  ga_map_grade_factor(ga, in, out, 2);
}

void ga_ep_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  int64_t left_term_index;
  int64_t right_term_index;
  ga_blade_id_t left_id;
  ga_blade_id_t right_id;
  ga_blade_id_t id;
  int64_t left_bit_index[ga_max_dimensions];
  int64_t right_bit_index[ga_max_dimensions];
  int64_t left_bit_index_count;
  int64_t right_bit_index_count;
  double sign;
  double coeff;

  ga_accumulator_clear(accumulator);

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    left_id = left->term[left_term_index].id;
    ga_id_to_bit_indices(ga, left_id, left_bit_index, &left_bit_index_count);

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      right_id = right->term[right_term_index].id;
      if (!(left_id & right_id)) {
        id = left_id | right_id;
        sign = 1.0;
        if (right_id) {
          ga_id_to_bit_indices(ga, right_id, right_bit_index, &right_bit_index_count);
          sign = ga_sign_sorted(left_bit_index, left_bit_index_count, right_bit_index, right_bit_index_count);
        }
        coeff = sign * left->term[left_term_index].coeff * right->term[right_term_index].coeff;
        ga_accumulator_add(accumulator, id, coeff);
      }
      right_term_index += 1;
    }

    left_term_index += 1;
  }

  ga_accumulator_to_mv(ga, accumulator, out);
}

void ga_gp_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  int64_t left_term_index;
  int64_t right_term_index;
  int64_t left_bit_index[ga_max_dimensions];
  int64_t right_bit_index[ga_max_dimensions];
  int64_t merged_bit_index[ga_max_dimensions * 2];
  int64_t left_bit_index_count;
  int64_t right_bit_index_count;
  int64_t merged_bit_index_count;
  int64_t index;
  int64_t changed;
  int64_t row_index;
  int64_t column_index;
  ga_blade_id_t id;
  double sign;
  double factor;
  double metric_value;
  double coeff;

  ga_accumulator_clear(accumulator);

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    ga_id_to_bit_indices(ga, left->term[left_term_index].id, left_bit_index, &left_bit_index_count);

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      ga_id_to_bit_indices(ga, right->term[right_term_index].id, right_bit_index, &right_bit_index_count);

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
        id = ga_id_from_bit_indices(merged_bit_index, merged_bit_index_count);
        ga_accumulator_add(accumulator, id, coeff);
      }

      right_term_index += 1;
    }

    left_term_index += 1;
  }

  ga_accumulator_to_mv(ga, accumulator, out);
}

void ga_ip_equal_grade(ga_t *ga, int64_t *left_bit_index, int64_t *right_bit_index, int64_t bit_index_count, ga_mv_t *out)
{
  double matrix[ga_max_dimensions][ga_max_dimensions];
  double scalar;
  int64_t row_index;
  int64_t column_index;
  int64_t sign_power;

  if (bit_index_count == 1) {
    ga_scalar(ga->metric[left_bit_index[0]][right_bit_index[0]], out);
    return;
  }

  row_index = 0;
  while (row_index < bit_index_count) {
    column_index = 0;
    while (column_index < bit_index_count) {
      matrix[row_index][column_index] = ga->metric[left_bit_index[row_index]][right_bit_index[column_index]];
      column_index += 1;
    }
    row_index += 1;
  }

  scalar = ga_determinant(matrix, bit_index_count);
  sign_power = (bit_index_count * (bit_index_count - 1)) / 2;
  if (sign_power & 1) {
    scalar = -scalar;
  }
  ga_scalar(-scalar, out);
}

void ga_choose_rec(int64_t *space, int64_t space_count, int64_t choose_count, int64_t space_index, int64_t chosen_count, int64_t *chosen, int64_t *combo_count, int64_t combo[ga_blade_count][ga_max_dimensions])
{
  int64_t copy_index;

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

double ga_sign_of_choice(int64_t *choice, int64_t choice_count)
{
  int64_t index;
  int64_t sum;

  sum = 0;
  index = 0;
  while (index < choice_count) {
    sum += choice[index] - index;
    index += 1;
  }
  return (sum & 1) ? -1.0 : 1.0;
}

void ga_ip_one(ga_t *ga, ga_mv_t *left, ga_mv_t *right, ga_accumulator_t *accumulator, ga_mv_t *out)
{
  int64_t left_term_index;
  int64_t right_term_index;
  int64_t left_bit_index[ga_max_dimensions];
  int64_t right_bit_index[ga_max_dimensions];
  int64_t left_bit_index_count;
  int64_t right_bit_index_count;
  int64_t space[ga_max_dimensions];
  int64_t chosen[ga_max_dimensions];
  int64_t combo[ga_blade_count][ga_max_dimensions];
  int64_t combo_count;
  int64_t combo_index;
  int64_t combo_position;
  int64_t remainder_bit_index[ga_max_dimensions];
  int64_t remainder_bit_index_count;
  int64_t index;
  int64_t found;
  ga_blade_id_t id;
  double sign;
  double coeff;
  ga_term_t scalar_storage[1];
  ga_mv_t scalar_mv;

  scalar_mv.term = scalar_storage;
  scalar_mv.term_count = 0;

  ga_accumulator_clear(accumulator);

  left_term_index = 0;
  while (left_term_index < left->term_count) {
    ga_id_to_bit_indices(ga, left->term[left_term_index].id, left_bit_index, &left_bit_index_count);

    right_term_index = 0;
    while (right_term_index < right->term_count) {
      ga_id_to_bit_indices(ga, right->term[right_term_index].id, right_bit_index, &right_bit_index_count);

      if (right_bit_index_count >= left_bit_index_count) {
        if (right_bit_index_count == left_bit_index_count) {
          ga_ip_equal_grade(ga, left_bit_index, right_bit_index, left_bit_index_count, &scalar_mv);
          if (scalar_mv.term_count) {
            coeff = left->term[left_term_index].coeff * right->term[right_term_index].coeff * scalar_mv.term[0].coeff;
            ga_accumulator_add(accumulator, 0, coeff);
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

            id = ga_id_from_bit_indices(remainder_bit_index, remainder_bit_index_count);
            coeff *= left->term[left_term_index].coeff * right->term[right_term_index].coeff * sign;
            ga_accumulator_add(accumulator, id, coeff);
            combo_index += 1;
          }
        }
      }

      right_term_index += 1;
    }

    left_term_index += 1;
  }

  ga_accumulator_to_mv(ga, accumulator, out);
}

status_t ga_inverse(ga_t *ga, ga_mv_t *in, ga_accumulator_t *accumulator, ga_mv_t *reverse_mv, ga_mv_t *denominator_mv, ga_mv_t *out)
{
  status_declare;
  int64_t term_index;
  double denominator;

  denominator = 0.0;

  ga_reverse(ga, in, reverse_mv);
  ga_gp_one(ga, in, reverse_mv, accumulator, denominator_mv);

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

  ga_scale(reverse_mv, 1.0 / denominator, out);

exit:
  status_return;
}

status_t ga_init_full(ga_t *ga, double metric[ga_max_dimensions][ga_max_dimensions], int64_t dimensions, int64_t is_conformal)
{
  status_declare;
  int64_t row_index;
  int64_t column_index;
  int64_t base_dimensions;

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
  } else {
    ga->null_vector_start = 0;
    ga->id_null = 0;
  }

  if (is_conformal) {
    ga->no_bit_index = dimensions - 2;
    ga->ni_bit_index = dimensions - 1;
    ga->no_id = (((ga_blade_id_t) 1)) << ga->no_bit_index;
    ga->ni_id = (((ga_blade_id_t) 1)) << ga->ni_bit_index;
    ga->no_index = ga->no_bit_index + 1;
    ga->ni_index = ga->ni_bit_index + 1;
  } else {
    ga->no_bit_index = 0;
    ga->ni_bit_index = 0;
    ga->no_id = 0;
    ga->ni_id = 0;
    ga->no_index = 0;
    ga->ni_index = 0;
  }

exit:
  if (status_is_failure) {
    ga_uninit(ga);
  }
  status_return;
}

status_t ga_init_diagonal(ga_t *ga, double *diagonal_metric, int64_t dimensions, int64_t is_conformal)
{
  status_declare;
  double full_metric[ga_max_dimensions][ga_max_dimensions];

  ga_flat_metric_to_full(diagonal_metric, dimensions, full_metric);
  status_require((ga_init_full(ga, full_metric, dimensions, is_conformal)));

exit:
  status_return;
}
