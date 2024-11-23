class sph_ga
  constructor: (metric, options = {}) ->
    @n = metric.length
    @is_conformal = !!options.conformal
    unless Array.isArray metric[0]
      if @is_conformal
        @n += 2
        metric.push 0, 0
      metric = @flat_metric_to_full metric, @n
      if @is_conformal
        metric[@n - 1][@n - 2] = -1
        metric[@n - 2][@n - 1] = -1
    [@is_symmetric, @is_diagonal, @null_vectors] = @metric_properties metric, @n
    unless @is_symmetric
      throw new Error "the metric must be symmetric in the non-null part"
    if @is_conformal
      unless 2 == @null_vectors
        throw new Error "only two null vectors are allowed with \"conformal: true\". use a custom metric instead"
      unless @is_diagonal
        throw new Error "only diagonal metrics are allowed with \"conformal: true\". use a custom metric instead"
    if @null_vectors
      @null_vector_start = @n - @null_vectors
      @id_null = @id_from_indices [(@null_vector_start + 1)..@n]
    @metric = metric
    @pseudoscalar_id = (1 << @n) - 1
    @ip_metric = @ip_metric_f()
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

  array_product: (a) -> a.reduce ((b, a) -> a * b), 1
  array_sum: (a) -> a.reduce ((b, a) -> a + b), 1
  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i, coeff) -> [@basis_blade(i, (if coeff? then coeff else 1))]
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a
  s: (coeff) -> [[0, coeff, 0]]
  id_from_indices: (indices) -> indices.reduce ((id, i) -> id |= 1 << (i - 1)), 0
  mv: (terms) -> @blade indices, coeff for [indices, coeff] in terms
  map_grade_factor: (a, f) -> ([id, coeff * f(grade), grade] for [id, coeff, grade] in a)
  involute: (a) -> @map_grade_factor a, (grade) -> (-1) ** grade
  scale: (mv, a) -> ([id, coeff * a, grade] for [id, coeff, grade] in mv)
  negate: (a) -> @scale a, -1
  conjugate: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade + 1)) >> 1)
  reverse: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade - 1)) >> 1)
  add: (a, b) -> @combine a, b, 1
  subtract: (a, b) -> @combine a, b, -1
  sp: (a, b) -> @gp @gp(a, b), @inverse(a)
  pseudoscalar: -> [@blade([1..@n], 1)]
  coeffs_add: (coeffs, id, coeff, grade) -> if coeffs[id]? then coeffs[id][0] += coeff else coeffs[id] = [coeff, grade]
  grade: (a) -> a[a.length - 1][2]
  blade_id: (a) -> a[0]
  blade_coeff: (a) -> a[1]
  blade_grade: (a) -> a[2]
  get: (a, id) -> for b in a then return b if id == b[0]
  bitcount_cache: {}
  id_bit_indices_cache: {}
  null_scalar: [[0, 0, 0]]

  blade: (indices, coeff) ->
    if indices[0] then [@id_from_indices(indices), coeff, indices.length]
    else [@id_from_indices(indices.slice(1)), coeff, indices.length - 1]

  coeffs_to_mv: (coeffs) ->
    a = ([parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0)
    if a.length then a else [[0, 0, 0]]

  id_bit_indices: (id) ->
    return @id_bit_indices_cache[id] if id of @id_bit_indices_cache
    a = (i for i in [0...@n] when id & (1 << i))
    @id_bit_indices_cache[id] = a
    a

  flat_metric_to_full: (metric, n) ->
    a = Array n
    for i in [0...n]
      b = Array(n).fill 0
      b[i] = metric[i]
      a[i] = b
    a

  bitcount: (a) ->
    return @bitcount_cache[a] if a of @bitcount_cache
    n = 0
    b = a
    while b != 0
      b &= b - 1
      n += 1
    @bitcount_cache[a] = n
    n

  count_inversions_sorted: (a, b) ->
    c = 0
    i = 0
    j = 0
    while i < a.length and j < b.length
      if a[i] <= b[j] then i += 1
      else
        c += a.length - i
        j += 1
    c

  determinant_generic = (a, n) ->
    return a[0][0] if n == 1
    b = 0
    for j in [0...n]
      c = (a[i][0...j].concat a[i][j + 1...]) for i in [1...n]
      sign = if 0 == j % 2 then 1 else -1
      b += sign * a[0][j] * determinant_generic(c)
    b

  determinant: (matrix) ->
    n = matrix.length
    if 1 == n then matrix[0][0]
    else if 2 == n then matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    else if 3 == n
      matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    else if 4 == n
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
    else @determinant_generic matrix, n

  metric_properties: (metric, n) ->
    null_vectors = 0; first_zero = n; is_diagonal = true
    for i in [0...n]
      di = metric[i][i]
      if di == 0
        null_vectors += 1; first_zero = i if first_zero == n
      else
        if first_zero < n or not di in [1, -1]
          is_diagonal = false; break
      for j in [(i + 1)...n]
        if i < first_zero
          if metric[i][j] != metric[j][i]
            return [false, false, null_vectors]
          if metric[i][j] != 0
            is_diagonal = false; break
      break unless is_diagonal
    for k in [first_zero...n]
      if metric[k][k] != 0
        return [false, false, null_vectors]
    [true, is_diagonal, null_vectors]

  indices_null_start: (a) ->
    b = a.findLastIndex (i) => i < @null_vector_start
    if -1 == b then 0 else b + 1

  indices_split_null: (a) ->
    start = @indices_null_start a
    [a.slice(0, start), a.slice(start)]

  blade_string: (a) ->
    [id, coeff, grade] = a
    if id then base = "e" + @id_bit_indices(id).map((b) -> b + 1).join "_"
    else base = "s"
    "#{coeff}#{base}"

  mv_string: (a) -> (@blade_string b for b in a).join " + "

  ip_metric_f: ->
    if @is_conformal
      (indices) ->
        switch indices.length - @indices_null_start(indices)
          when 2 then -1
          when 1 then 0
          else @array_product(@metric[i][i] for i in indices)
    else if @is_diagonal
      (indices) ->
        [left, right] = @indices_split_null indices
        if right.length > 0 then 0
        else
          if left.length > 0 then @array_product(@metric[i][i] for i in left) else 1
    else
      (indices) ->
        return 1 if !indices.length
        [left, right] = @indices_split_null indices
        if !right.length then @determinant((@metric[i][j] for j in left) for i in left)
        else if !left.length then 0
        else @determinant((@metric[i][j] for j in indices) for i in indices)

  ip: (a, b) ->
    coeffs = {}
    if a.length == 1 && !a[0][0]
      return if b.length == 1 && !b[0][0] then @null_scalar else ([id, coeff * a[0][1], grade] for [id, coeff, grade] in b)
    else if b.length == 1 && !b[0][0]
      return @null_scalar
    for [id_a, coeff_a, grade_a] in a
      indices_a = @id_bit_indices id_a
      for [id_b, coeff_b, grade_b] in b
        if id_a == id_b
          if metric_factor = @ip_metric indices_a
            @coeffs_add coeffs, 0, coeff_a * coeff_b * metric_factor, 0
        else if grade_a <= grade_b && id_a == (id_b & id_a)
          id_c = id_a ^ id_b
          indices_c = @id_bit_indices id_c
          if metric_factor = @ip_metric indices_c
            sign = (-1) ** @count_inversions_sorted indices_a, indices_c
            @coeffs_add coeffs, id_c, coeff_a * coeff_b * sign * metric_factor, grade_b - grade_a
        else if !((id_a | id_b) & ~@id_null)
          if metric_factor = @ip_metric @id_bit_indices id_a | id_b
            @coeffs_add coeffs, 0, coeff_a * coeff_b * metric_factor, 0
    @coeffs_to_mv coeffs

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      indices_a = @id_bit_indices id_a
      len_a = indices_a.length
      for [id_b, coeff_b, grade_b] in b when !(id_a & id_b)
        if !grade_a  # combinations involving scalars
          if grade_b then @coeffs_add coeffs, id_b, coeff_a * coeff_b, grade_b
          continue
        else if !grade_b
          @coeffs_add coeffs, id_a, coeff_a * coeff_b, grade_a
          continue
        id_c = id_a | id_b
        indices_b = @id_bit_indices id_b
        inversions = 0
        i = 0
        j = 0
        len_b = indices_b.length
        while i < len_a and j < len_b
          if indices_a[i] < indices_b[j] then i += 1
          else
            inversions += len_a - i
            j += 1
        coeff = (-1) ** inversions * coeff_a * coeff_b
        @coeffs_add coeffs, id_c, coeff, grade_a + grade_b
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
