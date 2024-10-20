class simple_ga
  constructor: (metric) ->
    @n = metric.length
    @metric = metric
    @blade_count = 2 ** @n

  blade_id: (grade) -> [1 << grade, grade]
  basis: (grade) -> [@blade_id(grade), [1]]

  grade: (a) ->
    # blade_id -> integer
    # determine grade by counting the number of set bits.
    n = 0
    while a != 0
      a &= a - 1
      n += 1
    n

  antisymmetric_sign: (a, b) ->
    # blade_id blade_id -> integer
    bits_a = a[0]
    bits_b = b[0]
    count = 0
    count_b = 0
    for k in [0...@n - 1]
      count += count_b if (bits_a >> k) & 1
      count_b += 1 if (bits_b >> k) & 1
    if (count & 1) == 0 then 1 else -1

  metric_sign: (a, b) ->
    # blade_id blade_id -> integer
    common_bases = a & b
    sign = 1
    i = 0
    while common_bases != 0
      sign *= @metric[i] if common_bases & 1
      common_bases >>= 1
      i += 1
    sign

  gp: (a, b) ->
    # multivector multivector -> multivector
    result = {}
    [a_ids, a_coeffs] = a
    [b_ids, b_coeffs] = b
    for i_a in [0...a_ids.length]
      i = a_ids[i_a]
      ai = a_coeffs[i_a]
      for i_b in [0...b_ids.length]
        j = b_ids[i_b]
        bj = b_coeffs[i_b]
        blade_id = i ^ j
        sign = @antisymmetric_sign(i, j) * @metric_sign(i, j)
        res_coeff = sign * ai * bj
        result[blade_id] = (result[blade_id] or 0) + res_coeff
    result_ids = []
    result_coeffs = []
    for blade_id, coeff of result when coeff != 0
      result_ids.push parseInt blade_id
      result_coeffs.push coeff
    [result_ids, result_coeffs]

  # Shared primitive function for grade-dependent sign
  apply_grade_sign: (a, sign_function) ->
    # multivector function -> multivector
    [ids, coeffs] = a
    new_ids = ids.slice(0)
    new_coeffs = []
    for i in [0...ids.length]
      blade = ids[i]
      coeff = coeffs[i]
      grade = @grade blade
      sign = sign_function(grade)
      new_coeffs[i] = coeff * sign
    [new_ids, new_coeffs]

  reverse: (a) ->
    # multivector -> multivector
    # reverse the order of basis vectors in each blade.
    # each blade's coefficient is multiplied by (-1 ** (k * (k - 1) / 2)), where k is the grade.
    sign_function = (grade) ->
      exponent = (grade * (grade - 1)) >> 1  # Integer division by 2
      (-1) ** exponent
    @apply_grade_sign a, sign_function

  involute: (a) ->
    # multivector -> multivector
    # changes the sign of blades based on their grade.
    # each blades coefficient is multiplied by (-1) ** k
    sign_function = (grade) -> (-1) ** grade
    @apply_grade_sign a, sign_function

  conjugate: (a) ->
    # multivector -> multivector
    # combines the reverse and involute operations.
    # each blades coefficient is multiplied by (-1) ** (k * (k + 1) / 2)
    sign_function = (grade) ->
      exponent = (grade * (grade + 1)) >> 1
      (-1) ** exponent
    @apply_grade_sign a, sign_function

  combine: (a, b, scalar = 1) ->
    # multivector multivector number -> multivector
    result = {}
    [a_ids, a_coeffs] = a
    [b_ids, b_coeffs] = b
    # add coefficients from a
    for i in [0...a_ids.length]
      blade_id = a_ids[i]
      coeff = a_coeffs[i]
      result[blade_id] = coeff
    # combine coefficients from b with scalar multiplier
    for i in [0...b_ids.length]
      blade_id = b_ids[i]
      coeff = scalar * b_coeffs[i]
      if blade_id of result then result[blade_id] += coeff
      else result[blade_id] = coeff
    result_ids = []
    result_coeffs = []
    for blade_id, coeff of result when coeff != 0
      result_ids.push parseInt blade_id
      result_coeffs.push coeff
    [result_ids, result_coeffs]

  add: (a, b) -> @combine a, b, 1
  subtract: (a, b) -> @combine a, b, -1

  ip: (a, b) ->
    # multivector multivector -> multivector
    result = {}
    [a_ids, a_coeffs] = a
    [b_ids, b_coeffs] = b
    for i_a in [0...a_ids.length]
      i = a_ids[i_a]
      ai = a_coeffs[i_a]
      bits_i = i[0]
      grade_i = i[1]
      for i_b in [0...b_ids.length]
        j = b_ids[i_b]
        bj = b_coeffs[i_b]
        bits_j = j[0]
        grade_j = j[1]
        if grade_i <= grade_j
          bits_blade_id = bits_i ^ bits_j
          grade_res = grade_j - grade_i
          grade_blade_id = grade_res
          blade_id = [bits_blade_id, grade_blade_id]
          sign = @antisymmetric_sign(i, j) * @metric_sign(i, j)
          res_coeff = sign * ai * bj
          key = bits_blade_id
          result[key] = [blade_id, 0] if key not in result
          result[key][1] += res_coeff
    result_ids = []
    result_coeffs = []
    for key, [blade_id, coeff] of result when coeff != 0
      result_ids.push blade_id
      result_coeffs.push coeff
    [result_ids, result_coeffs]

  ep: (a, b) ->
    # multivector multivector -> multivector
    result = {}
    [a_ids, a_coeffs] = a
    [b_ids, b_coeffs] = b
    for i_a in [0...a_ids.length]
      i = a_ids[i_a]
      ai = a_coeffs[i_a]
      grade_a = @grade i
      for i_b in [0...b_ids.length]
        j = b_ids[i_b]
        bj = b_coeffs[i_b]
        grade_b = @grade j
        blade_id = i ^ j
        grade_res = grade_a + grade_b
        if @grade(blade_id) == grade_res
          sign = @antisymmetric_sign(i, j) * @metric_sign(i, j)
          res_coeff = sign * ai * bj
          result[blade_id] = (result[blade_id] or 0) + res_coeff
    result_ids = []
    result_coeffs = []
    for blade_id, coeff of result when coeff != 0
      result_ids.push parseInt blade_id
      result_coeffs.push coeff
    [result_ids, result_coeffs]

  inverse: (a) ->
    # multivector -> multivector
    a_reverse = @reverse a
    denom_mv = @gp a, a_reverse
    denom = 0
    [ids, coeffs] = denom_mv
    for i in [0...ids.length]
      if ids[i] == 0  # scalar part has blade_id 0
        denom = coeffs[i]
        break
    if denom == 0
      throw new Error "multivector is not invertible (denominator is zero)."
    # compute inverse: a_inverse = a_reverse / denom
    a_inverse = [a_reverse[0], []]
    for i in [0...a_reverse[1].length]
      a_inverse[1][i] = a_reverse[1][i] / denom
    a_inverse

  sp: (a, b) ->
    # multivector multivector -> multivector
    # compute the sandwich product: a * b * a ** -1
    a_inverse = @inverse(a)
    temp = @gp a, b
    @gp temp, a_inverse

  pseudoscalar: ->
    # -> multivector
    ps = @basis 0
    ps = @ep ps, @basis(i) for i in [1...@n]
    ps
