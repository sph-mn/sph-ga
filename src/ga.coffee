class sph_ga
  constructor: (metric, options = {}) ->
    @n = metric.length
    @is_conformal = !!options.conformal
    if Array.isArray metric[0]
      [@is_orthonormal, @is_quasi_diagonal, @null_vectors] = @full_metric_properties metric, @n
    else
      [@is_orthonormal, @is_quasi_diagonal, @null_vectors] = @flat_metric_properties metric, @n
      metric = @flat_metric_to_full metric, @n
      if @is_conformal
        metric[@n - 2][@n - 1] = -1
        metric[@n - 1][@n - 2] = -1
    @metric = metric
    @pseudoscalar_id = (1 << @n) - 1
    if @is_conformal
      @eo_bit_index = @n - 2
      @ei_bit_index = @n - 1
      @eo_id = 1 << @eo_bit_index
      @ei_id = 1 << @ei_bit_index
      @eo_index = @eo_bit_index + 1
      @ei_index = @ei_bit_index + 1
      @eo = (coeff) -> [[@eo_id, coeff, 1]]
      @ei = (coeff) -> [[@ei_id, coeff, 1]]
      point = (euclidean_coeffs) ->
        ei_coeff = 0.5 * @array_sum(@array_product(euclidean_coeffs))
        @vector [0].concat(euclidean_coeffs).concat([1, ei_coeff])

  flat_metric_properties: (metric, n) ->
    null_vectors = 0
    is_orthonormal = true
    for a in metric
      if a then is_orthonormal = false unless a is 1 or a is -1
      else
        null_vectors += 1
        is_orthonormal = false
    [is_orthonormal, true, null_vectors]

  full_metric_properties: (metric, n) ->
    null_vectors = 0
    is_symmetric = true
    is_orthonormal = true
    is_quasi_diagonal = true
    for i in [0...n]  # check symmetry
      for j in [i + 1...n]
        unless metric[i][j] is metric[j][i]
          return [false, false, null_vectors]
    for i in [0...n]
      diagonal = metric[i][i]
      if !diagonal
        null_vectors += 1
        is_orthonormal = false
        is_quasi_diagonal = false
        continue
      unless diagonal is 1 or diagonal is -1
        is_orthonormal = false
      for j in [0...n]
        next if i is j
        unless metric[i][j] is 0
          is_quasi_diagonal = false
          is_orthonormal = false
          break
      break unless is_quasi_diagonal or is_orthonormal
    [is_orthonormal, is_quasi_diagonal, null_vectors]

  flat_metric_to_full: (metric, n) ->
    a = Array n
    for i in [0...n]
      b = Array(n).fill 0
      b[i] = metric[i]
      a[i] = b
    a

  array_product: (a) -> a.reduce ((b, a) -> a * b), 1
  array_sum: (a) -> a.reduce ((b, a) -> a + b), 1
  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i, coeff) -> [@basis_blade(i, (if coeff? then coeff else 1))]
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a
  s: (coeff) -> [[0, coeff, 0]]
  id_from_indices: (indices) -> indices.reduce ((id, i) -> id |= 1 << (i - 1)), 0
  blade: (indices, coeff) ->
    if indices[0] then [@id_from_indices(indices), coeff, indices.length]
    else [@id_from_indices(indices.slice(1)), coeff, indices.length - 1]
  mv: (terms) -> @blade indices, coeff for [indices, coeff] in terms
  map_grade_factor: (a, f) -> ([id, coeff * f(grade), grade] for [id, coeff, grade] in a)
  involute: (a) -> @map_grade_factor a, (grade) -> (-1) ** grade
  conjugate: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade + 1)) >> 1)
  reverse: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade - 1)) >> 1)
  add: (a, b) -> @combine a, b, 1
  subtract: (a, b) -> @combine a, b, -1
  sp: (a, b) -> @gp @gp(a, b), @inverse(a)
  pseudoscalar: -> [@blade([1..@n], 1)]
  coeffs_add: (coeffs, id, coeff, grade) -> if coeffs[id]? then coeffs[id][0] += coeff else coeffs[id] = [coeff, grade]
  coeffs_to_mv: (coeffs) -> [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0
  bitcount_cache: {}
  grade: (a) -> a[a.length - 1][2]
  blade_id: (a) -> a[0]
  blade_coeff: (a) -> a[1]
  blade_grade: (a) -> a[2]
  get: (a, id) -> for b in a then return b if id == b[0]

  bitcount: (a) ->
    return @bitcount_cache[a] if a of @bitcount_cache
    n = 0
    b = a
    while b != 0
      b &= b - 1
      n += 1
    @bitcount_cache[a] = n
    n

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

  ep_sign: (id_a, id_b) ->
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

  indices_cache: {}

  get_basis_indices: (id) ->
    return @indices_cache[id] if id of @indices_cache
    indices = []
    for i in [0...@n]
      if (id & (1 << i)) != 0
        indices.push i
    @indices_cache[id] = indices
    indices

  count_inversions_sorted: (a, b) ->
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

  ip_determinant_metric: (id, indices) -> @ip_metric(id, indices)

  ip_full_metric: (indices) ->
    if @is_conformal and ([ @eo_id, @ei_id ].some((id) -> indices.includes(id)))
      # special case: Null vectors in differing blade ids
      @ip_metric(@eo_id | @ei_id, indices)
    else @ip_metric 0, indices

  ip: (a, b) ->
    coeffs = {}
    if 1 == a.length  # scalar ⋅ b
      return ([id, coeff * a[0][1] * @ip_metric(id), grade] for [id, coeff, grade] in b) if !a[0][0]
    if 1 == b.length  # a ⋅ scalar
      return ([id, coeff * b[0][1] * @ip_metric(id), grade] for [id, coeff, grade] in a) if !b[0][0]
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
          sign = (-1) ** count_inversions_sorted indices_a, @get_basis_indices id_c
          coeff = coeff_a * coeff_b * sign * ip_full_metric indices_c
          @coeffs_add coeffs, id_a, coeff, grade_b - grade_a
    @coeffs_to_mv coeffs

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        continue if (id_a & id_b) != 0
        id_c = id_a | id_b
        coeff = coeff_a * coeff_b * @ep_sign id_a, idb
        @coeffs_add coeffs, id_c, coeff, grade_a + grade_b if coeff != 0
    @coeffs_to_mv coeffs

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

if typeof module isnt "undefined" and module.exports then module.exports = sph_ga
else window.sph_ga = sph_ga
