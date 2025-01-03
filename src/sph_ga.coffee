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
      null_vector_start = @n - @null_vectors
      @id_null = @id_from_indices [(null_vector_start + 1)..@n]
    @pseudoscalar_id = (1 << @n) - 1
    @metric = metric
    if @is_diagonal
      @ip_metric = (indices) -> @array_product(@metric[i][i] for i in indices)
    else
      @ip_metric = (indices) ->
        return 1 if !indices.length
        @determinant((@metric[i][j] for j in indices) for i in indices)
    if @is_conformal
      @no_bit_index = @n - 2
      @ni_bit_index = @n - 1
      @no_id = 1 << @no_bit_index
      @ni_id = 1 << @ni_bit_index
      @no_index = @no_bit_index + 1
      @ni_index = @ni_bit_index + 1
      @no = (coeff) -> [[@no_id, coeff, 1]]
      @ni = (coeff) -> [[@ni_id, coeff, 1]]
      rotation_axis_combinations = (n) ->
        combinations = []
        for i in [0...n]
          for j in [i+1...n]
            combinations.push([i + 1, j + 1])
        combinations
      @rotation_axes = rotation_axis_combinations @n - 2
      @rotor = (coeffs) ->
        [first, rest...] = coeffs
        @mv [[[0], first]].concat([@rotation_axes[i], a] for a, i in rest)
      @point = (euclidean_coeffs) ->
        ni_coeff = @array_sum(a * a for a in euclidean_coeffs) / 2
        @vector [0].concat(euclidean_coeffs).concat([1, ni_coeff])
      @normal = @vector [0].concat Array(@n).fill 1 / Math.sqrt @n

  add_one: (a, b) -> @combine a, b, 1
  array_diff: (a, b) -> a.filter (c) -> !(c in b)
  array_product: (a) -> a.reduce ((b, a) -> a * b), 1
  array_sum: (a) -> a.reduce ((b, a) -> a + b), 0
  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i, coeff=1) -> [@basis_blade(i, coeff)]
  blade_coeff: (a) -> a[1]
  blade_grade: (a) -> a[2]
  blade_id: (a) -> a[0]
  coeffs_add: (coeffs, id, coeff, grade) -> if coeffs[id]? then coeffs[id][0] += coeff else coeffs[id] = [coeff, grade]
  conjugate: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade + 1)) >> 1)
  get: (a, id) -> for b in a then return b if id == b[0]
  grade: (a) -> a[a.length - 1][2]
  id_bit_indices_cache: {}
  id_from_bit_indices: (indices) -> indices.reduce ((id, i) -> id |= 1 << i), 0
  id_from_indices: (indices) -> indices.reduce ((id, i) -> id |= 1 << (i - 1)), 0
  id_grade_cache: {}
  involute: (a) -> @map_grade_factor a, (grade) -> (-1) ** grade
  map_grade_factor: (a, f) -> ([id, coeff * f(grade), grade] for [id, coeff, grade] in a)
  mv: (terms) -> @blade indices, coeff for [indices, coeff] in terms
  negate: (a) -> @scale a, -1
  null_scalar: [[0, 0, 0]]
  pseudoscalar: -> [@blade([1..@n], 1)]
  reverse: (a) -> @map_grade_factor a, (grade) -> (-1) ** ((grade * (grade - 1)) >> 1)
  scale: (mv, a) -> ([id, coeff * a, grade] for [id, coeff, grade] in mv)
  s: (coeff) -> [[0, coeff, 0]]
  sp: (a, b) -> @gp @gp(a, b), @inverse(a)
  subtract_one: (a, b) -> @combine a, b, -1
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a

  blade: (indices, coeff) ->
    if indices[0] then [@id_from_indices(indices), coeff, indices.length]
    else [@id_from_indices(indices.slice(1)), coeff, indices.length - 1]

  coeffs_to_mv: (coeffs) ->
    a = ([parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0)
    if a.length then a else [[0, 0, 0]]

  id_indices: (id) -> if id then (1 + a for a in @id_bit_indices(id)) else [0]

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

  id_grade: (a) ->
    return @id_grade_cache[a] if a of @id_grade_cache
    n = 0
    b = a
    while b != 0
      b &= b - 1
      n += 1
    @id_grade_cache[a] = n
    n

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

  for_each_combination: (array, n, f) ->
    generate = (prefix, rest, k, i) ->
      if k == 0
        f prefix, i
        return i + 1
      return i if rest.length == 0
      [first, rest...] = rest
      i = generate prefix.concat([first]), rest, k - 1, i
      generate prefix, rest, k, i
    generate [], array, n, 0

  sign: (indices) ->
    c = 0
    for i in [0...indices.length]
      for j in [(i + 1)...indices.length]
        c += 1 if indices[i] > indices[j]
    (-1) ** c

  sign_sorted: (a, b) ->
    c = 0
    i = 0
    j = 0
    while i < a.length and j < b.length
      if a[i] <= b[j]
        i += 1
      else
        j += 1
    (-1) ** c

  merge_indices: (a) ->
    factor = 1
    changed = true
    while changed
      changed = false
      i = 0
      while i < a.length
        j = i + 1
        while j < a.length
          if a[i] != a[j]
            m = @metric[a[i]][a[j]]
            if m != 0
              factor *= m
              a.splice j, 1
              a.splice i, 1
              changed = true
              break
          j += 1
        break if changed
        i += 1 unless changed
    changed = true
    while changed
      changed = false
      i = 0
      while i < a.length
        j = i + 1
        while j < a.length
          if a[i] is a[j]
            m = @metric[a[i]][a[i]]
            factor *= m
            a.splice j, 1
            a.splice i, 1
            changed = true
            break
          j += 1
        break if changed
        i += 1 unless changed
    [a, factor]

  ip_one: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        continue if grade_b < grade_a  # Left contraction requires grade_b >= grade_a
        indices = @id_bit_indices(id_a).concat(@id_bit_indices(id_b)).sort()
        [indices, factor] = @merge_indices indices
        coeff = coeff_a * coeff_b * factor
        grade = indices.length
        continue if coeff == 0
        target_grade = grade_b - grade_a
        continue unless grade is target_grade and coeff != 0
        @coeffs_add coeffs, @id_from_bit_indices(indices), coeff, grade
    @coeffs_to_mv coeffs

  gp_one: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      indices_a = @id_bit_indices id_a
      for [id_b, coeff_b, grade_b] in b
        indices_ab = indices_a.concat @id_bit_indices id_b
        indices = indices_ab.slice().sort()
        sign = @sign indices_ab
        coeff = coeff_a * coeff_b
        [indices, factor] = @merge_indices indices
        coeff *= factor * sign
        continue unless coeff
        @coeffs_add coeffs, @id_from_bit_indices(indices), coeff, indices.length
    @coeffs_to_mv coeffs

  ep_one: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      indices_a = @id_bit_indices id_a
      for [id_b, coeff_b, grade_b] in b when !(id_a & id_b)
        id = id_a | id_b
        continue unless id
        sign = if id_b then @sign_sorted indices_a, @id_bit_indices id_b else 1
        @coeffs_add coeffs, id, sign * coeff_a * coeff_b, grade_a + grade_b
    @coeffs_to_mv coeffs

  ep: (...a) -> a.reduce (c, b) => @ep_one c, b
  ip: (...a) -> a.reduce (c, b) => @ip_one c, b
  gp: (...a) -> a.reduce (c, b) => @gp_one c, b
  add: (...a) -> a.reduce (c, b) => @add_one c, b
  subtract: (...a) -> a.reduce (c, b) => @subtract_one c, b

  combine: (a, b, scalar = 1) ->
    coeffs = {}
    coeffs[id] = [coeff, grade] for [id, coeff, grade] in a
    for [id, coeff, grade] in b
      if coeffs[id]? then coeffs[id][0] += coeff * scalar
      else coeffs[id] = [coeff * scalar, grade]
    c = ([parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0)
    if c.length then c else @null_scalar

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

  blade_to_string: (a) ->
    [id, coeff, grade] = a
    if id
      base = "e" + @id_bit_indices(id).map((b) -> b + 1).join "_"
      if 1 == coeff then base else coeff + base
    else coeff

  mv_to_string: (a) -> (@blade_to_string b for b in a).join " + "

  blade_from_string: (a) ->
    match = a.match /^(?:(\d+(?:\.\d+)?))?([a-z]+)?(?:([\d_]+))?$/
    return null unless match
    left_number = if match[1]? then parseFloat match[1] else null
    letters = if match[2]? then match[2] else null
    right_numbers = if match[3]? then (parseInt n for n in match[3].split "_") else null
    result = []
    coeff = if left_number? then left_number else 1
    indices = if right_numbers? then right_numbers else []
    if letters?
      switch letters
        when "no" then id = c3.no_id
        when "ni" then id = c3.ni_id
        else id = @id_from_indices right_numbers
    else id = 0
    [id, coeff, @id_grade id]

  mv_from_string: (a) -> @blade_from_string b for b in a.split " + "

if typeof module isnt "undefined" and module.exports then module.exports = sph_ga
else window.sph_ga = sph_ga
