class sph_ga
  constructor: (euclidean_metric, is_conformal=false) ->
    @n = euclidean_metric.length
    @metric = euclidean_metric
    @is_conformal = !!is_conformal
    if is_conformal
      @n += 2
      @metric.push 0, 0
      @eo_index = @n - 1
      @ei_index = @n
      @eo_bit_index = @eo_index - 1
      @ei_bit_index = @ei_index - 1
      @eo_id = 1 << @eo_bit_index
      @ei_id = 1 << @ei_bit_index
      @eo = (coeff) -> [[@eo_id, coeff, 1]]
      @ei = (coeff) -> [[@ei_id, coeff, 1]]
    @pseudoscalar_id = (1 << @n) - 1
    @metric_map = {}
    for i in [0...@n]
      @metric_map[1 << i] = @metric[i]

  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i, coeff) -> [@basis_blade(i, (if coeff? then coeff else 1))]
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a
  s: (coeff) -> [[0, coeff, 0]]
  blade: (indices, coeff) ->
    id = indices.reduce ((id, i) -> if !i then id else id |= 1 << (i - 1)), 0
    [id, coeff, @grade(id)]
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

  coeffs_mv: (coeffs) -> [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

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

  # Memoization cache for index retrieval
  indices_cache: {}

  get_indices: (id) ->
    return @indices_cache[id] if id of @indices_cache
    indices = []
    for i in [0...@n]
      if (id & (1 << i)) != 0
        indices.push i
    @indices_cache[id] = indices
    indices

  compute_ip_sign_and_metric: (id_a, id_b, shared) ->
    indices_a = @get_indices(id_a)
    indices_b = @get_indices(id_b)
    indices_shared = @get_indices(shared)
    # Positions of shared indices in b
    positions_in_b = [ indices_b.indexOf(idx) for idx in indices_shared ]
    swap_count = positions_in_b.reduce ((sum, pos) -> sum + pos), 0
    # Compute sign based on swap_count
    sign = (-1) ** swap_count
    # Compute metric value
    metric_val = 1
    for idx in indices_shared
      if @is_conformal and (idx == @eo_bit_index or idx == @ei_bit_index)
        metric_val = 0
        break
      else
        sign *= -1 if @metric[idx] < 0
        metric_val *= @metric_map[1 << idx]
    return { sign, metric_val }

  compute_ep_sign: (id_a, id_b) ->
    vectors_a = @get_indices(id_a)
    vectors_b = @get_indices(id_b)
    concatenated = vectors_a.concat(vectors_b)
    # Count the number of inversions (swaps needed to sort the list)
    count = 0
    for i in [0...concatenated.length]
      for j in [i + 1...concatenated.length]
        if concatenated[i] > concatenated[j]
          count += 1
    sign = if (count % 2) == 0 then 1 else -1
    sign

  ip_conformal_special_cases: (id_a, id_b, coeff_a, coeff_b, grade_b) ->
    return null unless @is_conformal
    # Inner product of eo and ei yields scalar -1
    if (id_a == @eo_id and id_b == @ei_id) or (id_a == @ei_id and id_b == @eo_id)
      coeff = -1 * coeff_a * coeff_b
      return { id_result: 0, coeff, grade_result: 0 }
    # Inner product of eo and (A ∧ ei) yields -A
    if id_a == @eo_id and (id_b & @ei_id) != 0 and grade_b > 1
      id_result = id_b ^ @ei_id
      grade_A = grade_b - 1
      sign = (-1) ** grade_A
      coeff = sign * coeff_a * coeff_b
      return { id_result, coeff, grade_result: grade_A }
    # Inner product of ei and (A ∧ eo) yields A
    if id_a == @ei_id and (id_b & @eo_id) != 0 and grade_b > 1
      id_result = id_b ^ @eo_id
      coeff = coeff_a * coeff_b  # Sign is positive
      grade_result = grade_b - 1
      return { id_result, coeff, grade_result }
    # Inner product of eo with eo or ei with ei is zero
    if (id_a == @eo_id and id_b == @eo_id) or (id_a == @ei_id and id_b == @ei_id)
      return { id_result: 0, coeff: 0, grade_result: 0 }
    return null

  ip: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        # Skip inner product if either operand is scalar
        continue if grade_a == 0 or grade_b == 0
        # Handle special cases
        special_case = @ip_conformal_special_cases(id_a, id_b, coeff_a, coeff_b, grade_b)
        if special_case?
          { id_result, coeff, grade_result } = special_case
          @coeffs_add(coeffs, id_result, coeff, grade_result) if coeff != 0
          continue
        # General case: Inner product is defined when grade_a <= grade_b
        continue if grade_a > grade_b
        grade_result = grade_b - grade_a
        # Compute shared basis vectors
        shared = id_a & id_b
        k = @grade(shared)
        # The grade of the geometric product
        grade_gp = grade_a + grade_b - 2 * k
        # Inner product contributes when grade_gp == grade_result
        continue if grade_gp != grade_result
        # Compute the resulting blade ID
        id_result = id_a ^ id_b
        # Compute sign and metric value
        { sign, metric_val } = @compute_ip_sign_and_metric(id_a, id_b, shared)
        continue if metric_val == 0  # Skip if metric value is zero
        # Compute the coefficient
        coeff = sign * metric_val * coeff_a * coeff_b
        @coeffs_add(coeffs, id_result, coeff, grade_result) if coeff != 0
    # If coeffs is empty, result is zero scalar
    @coeffs_add(coeffs, 0, 0, 0) if Object.keys(coeffs).length == 0
    @coeffs_mv(coeffs)

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        if (id_a & id_b) != 0 then continue  # Overlapping basis vectors result in zero
        id_result = id_a | id_b
        sign = @compute_ep_sign(id_a, id_b)
        coeff = sign * coeff_a * coeff_b
        @coeffs_add(coeffs, id_result, coeff, grade_a + grade_b) if coeff != 0
    @coeffs_mv(coeffs)

  gp: (a, b) ->
    # Compute the inner product
    ip_result = @ip(a, b)
    # Compute the exterior product
    ep_result = @ep(a, b)
    console.log "ip_result", ip_result
    console.log "ep_result", ep_result
    # Initialize a map to accumulate coefficients by blade ID
    coeffs = {}
    # Helper function to add terms to the coeffs map
    add_terms = (result) =>
      for [id, coeff, grade] in result
        @coeffs_add(coeffs, id, coeff, grade) if coeff != 0
    # Add inner product results
    add_terms(ip_result)
    # Add exterior product results
    add_terms(ep_result)
    # If no terms were added, return zero scalar
    if Object.keys(coeffs).length == 0
      return @coeffs_mv([[0, 0, 0]])
    # Otherwise, return the accumulated multivector
    return @coeffs_mv(coeffs)

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
