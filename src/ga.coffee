class sph_ga
  constructor: (metric, is_conformal, options) ->
    @n = metric.length
    @metric = if Array.isArray(metric[0]) then metric else @diagonal_metric metric
    @pseudoscalar_id = (1 << @n) - 1
    null_vectors = option.null_vector_indices || []
    @is_conformal = is_conformal
    if @is_conformal
      @eo_bit_index = @n - 2
      @ei_bit_index = @n - 1
      @eo_id = 1 << @eo_bit_index
      @ei_id = 1 << @ei_bit_index
      @eo_index = @eo_bit_index + 1
      @ei_index = @ei_bit_index + 1
      @eo = (coeff) -> [[@eo_id, coeff, 1]]
      @ei = (coeff) -> [[@ei_id, coeff, 1]]
      null_vector_indices.push @eo_index, @ei_index
      null_vector_indices.sort()
    @null_vectors = new Set null_vector_indices

  array_product = (a) -> a.reduce ((b, a) -> a * b), 1
  array_sum = (a) -> a.reduce ((b, a) -> a + b), 1

  diagonal_metric: (a) ->
    b = []
    for i in [0...@n]
      b[i] = []
      for j in [0...@n]
        b[i][j] = if i == j then a[i] else 0

  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i, coeff) -> [@basis_blade(i, (if coeff? then coeff else 1))]
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a
  s: (coeff) -> [[0, coeff, 0]]
  id_from_indices: (indices) -> indices.reduce ((id, i) -> id |= 1 << (i - 1)), 0
  blade: (indices, coeff) ->
    if indices[0] then [id_from_indices(indices), coeff, indices.length]
    else [id_from_indices(indices.slice(1)), coeff, indices.length - 1]
  mv: (terms) -> @blade indices, coeff for [indices, coeff] in terms
  apply_grade_sign: (a, f) -> ([id, coeff * f(grade), grade] for [id, coeff, grade] in a)
  involute: (a) -> @apply_grade_sign a, (grade) -> (-1) ** grade
  conjugate: (a) -> @apply_grade_sign a, (grade) -> (-1) ** ((grade * (grade + 1)) >> 1)
  reverse: (a) -> @apply_grade_sign a, (grade) -> (-1) ** ((grade * (grade - 1)) >> 1)
  add: (a, b) -> @combine a, b, 1
  subtract: (a, b) -> @combine a, b, -1
  sp: (a, b) -> @gp @gp(a, b), @inverse(a)
  pseudoscalar: -> [@blade([1..@n], 1)]

  coeffs_add: (coeffs, id, coeff, grade) ->
    if coeffs[id]? then coeffs[id][0] += coeff
    else coeffs[id] = [coeff, grade]

  coeffs_to_mv: (coeffs) -> [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  # Memoization cache for grade computations
  grade_cache: {}

  grade: (a) ->
    return @grade_cache[a] if a of @grade_cache
    n = 0
    b = a
    while b != 0
      b &= b - 1
      n += 1
    @grade_cache[a] = n
    n


  mv_blade: (a, id) ->
    for b in a
      return b if id == b[0]

  blade_id: (a) -> a[0]
  blade_coeff: (a) -> a[1]
  blade_grade: (a) -> a[2]

  determinant: (matrix) ->
    determinant_generic = (matrix) ->
      # Generic solution for larger matrices using Laplace expansion
      n = matrix.length
      return matrix[0][0] if n is 1
      det = 0
      for j in [0...n]
        submatrix = (matrix[i][0...j].concat matrix[i][j + 1...]) for i in [1...n]
        sign = if j % 2 is 0 then 1 else -1
        det += sign * matrix[0][j] * determinant_generic(submatrix)
      det
    n = matrix.length
    # Hardcoded formulas for small matrices
    if n is 1 then matrix[0][0]
    else if n is 2 then matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    else if n is 3
      matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    else if n is 4
      matrix[0][0] * (
        matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
        matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
        matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])) -
      matrix[0][1] * (
        matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
        matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
        matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])) +
      matrix[0][2] * (
        matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
        matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
        matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])) -
      matrix[0][3] * (
        matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
        matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
        matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]))
    else determinant_generic(matrix)

  compute_ep_sign: (id_a, id_b) ->
    vectors_a = @get_basis_indices(id_a)
    vectors_b = @get_basis_indices(id_b)
    concatenated = vectors_a.concat(vectors_b)
    # Count the number of inversions (swaps needed to sort the list)
    count = 0
    for i in [0...concatenated.length]
      for j in [i + 1...concatenated.length]
        if concatenated[i] > concatenated[j]
          count += 1
    sign = if (count % 2) == 0 then 1 else -1
    sign

  # Memoization cache for index retrieval
  indices_cache: {}

  get_basis_indices: (id) ->
    return @indices_cache[id] if id of @indices_cache
    indices = []
    for i in [0...@n]
      if (id & (1 << i)) != 0
        indices.push i
    @indices_cache[id] = indices
    indices

  count_inversions_sorted = (a, b) ->
    inversions = 0
    i = 0
    j = 0
    while i < a.length and j < b.length
      if a[i] <= b[j]
        i += 1
      else
        inversions += a.length - i
        j += 1
    inversions

  ip_metric = (id) ->
    indices = @get_basis_indices id # Extract basis indices for the blade
    a = 1
    for i in indices
      for j in indices
        a *= @metric[i][j] || 0 # Direct access into the metric array
    a

  metric_contribution = (id, indices) ->
    if @is_orthonormal
      if @is_conformal and (id & @eo_id or id & @ei_id)
        # Special case for conformal null vectors
        if id & @eo_id and id & @ei_id
          array_product([@metric[i][i] for i in indices]) * -1
        else 0
      else
        array_product([@metric[i][i] for i in indices])
    else
      # Non-orthonormal metric calculation
      determinant = 1
      for i in indices
        for j in indices
          determinant *= @metric[i][j]
      determinant

  ip_determinant_metric = (id, indices) ->
    metric_contribution(id, indices)

  ip_full_metric = (indices) ->
    if @is_conformal and ([ @eo_id, @ei_id ].some((id) -> indices.includes(id)))
      # Special case: Null vectors in differing blade IDs
      metric_contribution(@eo_id | @ei_id, indices)  # Use logical OR for shared context
    else
      metric_contribution 0, indices

  ip: (a, b) ->
    coeffs = {}
    if 1 == a.length  # scalar ⋅ b
      return ([id, coeff * a[0][1] * ip_metric(id), grade] for [id, coeff, grade] in b) if !a[0][0]
    if 1 == b.length  # a ⋅ scalar
      return ([id, coeff * b[0][1] * ip_metric(id), grade] for [id, coeff, grade] in a) if !b[0][0]
    for [id_a, coeff_a, grade_a] in a
      indices_a = @get_basis_indices id_a
      for [id_b, coeff_b, grade_b] in b when grade_a <= grade_b && (id_a & id_b) == id_a
        if id_b == id_a
          if 1 == grade_b && 1 == grade_a  # e_i ⋅ e_i
            i = indices_a[0]
            @coeffs_add coeffs, id_a, @metric[i][i], 1
          else  # a ⋅ a
            @coeffs_add coeffs, id_a, coeff_a * coeff_b * ip_determinant_metric(id_a, indices_a), grade_a
        else  # a ⋅ b
          id_c = id_b & ~id_a
          sign = (-1) ** count_inversions_sorted indices_a, @get_indices id_c
          coeff = coeff_a * coeff_b * sign * ip_full_metric indices_c
          @coeffs_add coeffs, id_a, coeff, grade_b - grade_a
    @coeffs_to_mv coeffs

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        if (id_a & id_b) != 0 then continue  # Overlapping basis vectors result in zero
        id_result = id_a | id_b
        sign = @compute_ep_sign(id_a, id_b)
        coeff = sign * coeff_a * coeff_b
        @coeffs_add(coeffs, id_result, coeff, grade_a + grade_b) if coeff != 0
    @coeffs_to_mv(coeffs)

  gp: (a, b) -> @s 0

  combine: (a, b, scalar = 1) ->
    coeffs = {}
    coeffs[id] = [coeff, grade] for [id, coeff, grade] in a
    for [id, coeff, grade] in b
      if coeffs[id]? then coeffs[id][0] += coeff * scalar
      else coeffs[id] = [coeff * scalar, grade]
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  inverse: (a) ->
    a_reverse = @reverse a
    denom_mv = @gp a, a_reverse
    denom = 0
    for [id, coeff, grade] in denom_mv
      if id == 0
        denom = coeff
        break
    if denom == 0 then throw new Error "multivector is not invertible (denominator is zero)."
    [parseInt(id), coeff / denom, grade] for [id, coeff, grade] in a_reverse

  point: (euclidean_coeffs) ->
    ei_coeff = 0.5 * (a * a for a in euclidean_coeffs).reduce(((b, a) -> b + a), 0)
    @vector [1].concat(euclidean_coeffs).concat([1, ei_coeff])

if typeof module isnt "undefined" and module.exports then module.exports = sph_ga
else window.sph_ga = sph_ga
