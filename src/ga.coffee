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
      @eo_id = 1 << (@eo_index - 1)
      @ei_id = 1 << (@ei_index - 1)
      @eo = (coeff) -> [[@eo_id, coeff, 1]]
      @ei = (coeff) -> [[@ei_id, coeff, 1]]
    @signature = @metric_to_signature @metric
    @pseudoscalar_id = (1 << @n) - 1

  basis_blade: (i, coeff) -> if i then [1 << (i - 1), coeff, 1] else [0, coeff, 0]
  basis: (i) -> [@basis_blade i, 1]
  vector: (coeffs) -> @basis_blade i, a for a, i in coeffs when a
  s: (coeff) -> [[0, coeff, 0]]
  blade: (indices, coeff) -> [indices.reduce(((id, i) -> if !i then id else id |= 1 << (i - 1)), 0), coeff, indices.length]
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

  metric_to_signature: (a) ->
    p = q = r = 0
    for i in a
      if i > 0 then p += 1
      else if i < 0 then q += 1
      else r += 1
    [p, q, r]

  grade: (a) ->
    n = 0
    while a != 0
      a &= a - 1
      n += 1
    n

  get_indices: (id) ->
    indices = []
    for i in [0...@n]
      if (id & (1 << i)) != 0
        indices.push i
    indices

  antisymmetric_sign: (id_a, id_b) ->
    vectors_a = @get_indices(id_a)
    vectors_b = @get_indices(id_b)
    # Concatenate vectors from both blades, vectors_b first
    concatenated = vectors_b.concat(vectors_a)
    # Count the number of inversions (swaps needed to sort the list)
    count = 0
    for i in [0...concatenated.length]
      for j in [i + 1...concatenated.length]
        if concatenated[i] > concatenated[j]
          count += 1
    sign = if (count % 2) == 0 then 1 else -1
    sign

  contraction_sign: (id_a, id_b) ->
    vectors_a = @get_indices(id_a)
    vectors_b = @get_indices(id_b)
    # Identify shared indices
    shared_indices = vectors_a.filter (idx) -> idx in vectors_b
    # Remove shared indices
    vectors_a_remain = vectors_a.filter (idx) -> not idx in shared_indices
    vectors_b_remain = vectors_b.filter (idx) -> not idx in shared_indices
    # Concatenate remaining vectors: vectors_a_remain first
    concatenated = vectors_a_remain.concat(vectors_b_remain)
    # Count inversions to sort concatenated list
    count = 0
    for i in [0...concatenated.length]
      for j in [i + 1...concatenated.length]
        if concatenated[i] > concatenated[j]
          count += 1
    sign = if (count % 2) == 0 then 1 else -1
    sign

  compute_metric: (id_a) ->
    metric_val = 1
    for i in [0...@n]
      if (id_a & (1 << i)) != 0
        # Exclude null vectors eo and ei
        if !(@is_conformal and (i == @eo_id or i == @ei_id))
          metric_val *= @metric[i]
    metric_val

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        # Handle scalar blades (grade 0)
        if grade_a == 0 and grade_b == 0
          # Scalar ^ Scalar = 0
          continue
        if grade_a == 0
          # Scalar ^ Blade = Blade scaled by scalar
          id = id_b
          grade = grade_b
          coeff = coeff_a * coeff_b
          if coeff != 0
            if coeffs[id]? then coeffs[id][0] += coeff
            else coeffs[id] = [coeff, grade]
          continue
        if grade_b == 0
          # Blade ^ Scalar = Blade scaled by scalar
          id = id_a
          grade = grade_a
          coeff = coeff_a * coeff_b
          if coeff != 0
            if coeffs[id]? then coeffs[id][0] += coeff
            else coeffs[id] = [coeff, grade]
          continue
        # For higher-grade blades, ensure no overlapping basis vectors
        if (id_a & id_b) != 0 then continue  # Overlapping basis vectors result in zero
        id = id_a | id_b  # Use bitwise OR to combine indices correctly
        sign = @antisymmetric_sign(id_a, id_b)
        coeff = sign * coeff_a * coeff_b
        if coeff != 0
          if coeffs[id]? then coeffs[id][0] += coeff
          else
            coeffs[id] = [coeff, grade_a + grade_b]
    # Convert to multivector format, filtering out zero coefficients
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  ip: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        # Skip inner product if either operand is scalar
        if grade_a == 0 or grade_b == 0 then continue

        # Inner product is defined when grade_a <= grade_b
        if grade_a > grade_b then continue

        grade_diff = grade_b - grade_a

        # Special case: Inner product of eo and ei should yield scalar -1
        if (@eo_id == id_a and @ei_id == id_b) or (@ei_id == id_a and @eo_id == id_b)
          coeff = -1 * coeff_a * coeff_b
          grade_result = 0  # Scalar
          @coeffs_add(coeffs, 0, coeff, grade_result, 0) if coeff != 0
          continue

        # Special case: Inner product of eo and (A ∧ ei) yields A
        if id_a == @eo_id and (id_b & @ei_id) != 0 and grade_b > 1
          # Remove ei from id_b to get id_result representing A
          id_result = id_b ^ @ei_id
          # Compute the sign due to swapping ei past A
          num_swaps = grade_b - 1
          sign = (-1) ** num_swaps
          # Compute the coefficient
          coeff = sign * coeff_a * coeff_b
          # The grade of the resulting blade
          grade_result = grade_b - 1
          @coeffs_add(coeffs, id_result, coeff, grade_result, 0) if coeff != 0
          continue

        # Special case: Inner product of ei and (A ∧ eo) yields A
        if id_a == @ei_id and (id_b & @eo_id) != 0 and grade_b > 1
          # Remove eo from id_b to get id_result representing A
          id_result = id_b ^ @eo_id
          # Compute the sign due to swapping eo past A
          num_swaps = grade_b - 1
          sign = (-1) ** num_swaps
          # Compute the coefficient
          coeff = sign * coeff_a * coeff_b
          # The grade of the resulting blade
          grade_result = grade_b - 1
          @coeffs_add(coeffs, id_result, coeff, grade_result, 0) if coeff != 0
          continue

        # General inner product for blades where grade_a <= grade_b
        # Check if blade A is a subset of blade B
        if (id_a & id_b) == id_a
          # Compute the resulting blade ID by removing blade A's bits from blade B
          id_result = id_b ^ id_a

          # Compute the sign based on the ordering of basis vectors
          sign = @antisymmetric_sign(id_a, id_b)

          # Compute metric value; assuming Euclidean for simplicity
          # Adjust if handling different metrics
          metric_val = 1
          for i in [0...@n]
            if (id_a & (1 << i)) != 0 and !(@is_conformal and (@eo_id == i or @ei_id == i))
              metric_val *= @metric[i]

          # Compute the coefficient
          coeff = sign * metric_val * coeff_a * coeff_b

          # The grade of the resulting blade
          grade_result = grade_diff

          @coeffs_add(coeffs, id_result, coeff, grade_result, 0) if coeff !=0

        # Enhancement 1: Contraction Logic for Partial Overlaps
        else
          # Identify shared basis vectors
          shared = id_a & id_b
          if shared !=0
            # Remove shared basis vectors from both blades
            id_a_remaining = id_a ^ shared
            id_b_remaining = id_b ^ shared

            # Compute the resulting blade ID by combining remaining parts
            id_result = id_a_remaining ^ id_b_remaining

            # Compute the sign based on the ordering of basis vectors during contraction
            sign = @antisymmetric_sign(shared, id_a_remaining)

            # Compute metric value; assuming Euclidean for simplicity
            metric_val =1
            for i in [0...@n]
              if (shared & (1 << i)) !=0 and !(@is_conformal and (@eo_id == i or @ei_id == i))
                metric_val *= @metric[i]

            # Compute the coefficient
            coeff = sign * metric_val * coeff_a * coeff_b

            # The grade of the resulting blade
            grade_result = grade_diff

            @coeffs_add(coeffs, id_result, coeff, grade_result, 0) if coeff !=0

        continue # finish the inner loop

    @coeffs_mv(coeffs)

  gp: (a, b) ->
    inner = @ip a, b
    outer = @ep a, b
    combined = {}
    for [id, coeff, grade] in inner
      @coeffs_add combined, id, coeff, grade
    for [id, coeff, grade] in outer
      @coeffs_add combined, id, coeff, grade
    @coeffs_mv combined

  combine: (a, b, scalar = 1) ->
    ##console.log "combine input:", a, b
    coeffs = {}
    coeffs[id] = [coeff, grade] for [id, coeff, grade] in a
    for [id, coeff, grade] in b
      if coeffs[id]? then coeffs[id][0] += coeff * scalar
      else coeffs[id] = [coeff * scalar, grade]
    ##console.log "combine output:", ([parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0)
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  inverse: (a) ->
    a_reverse = @reverse a
    console.log "inverse", a
    console.log "inverse reverse", a_reverse
    denom_mv = @gp a, a_reverse
    denom = 0
    for [id, coeff, grade] in denom_mv
      if id == 0
        denom = coeff
        break
    if denom == 0 then throw new Error "multivector is not invertible (denominator is zero)."
    #console.log ([parseInt(id), coeff / denom, grade] for [id, coeff, grade] in a_reverse)
    [parseInt(id), coeff / denom, grade] for [id, coeff, grade] in a_reverse

  point: (euclidean_coeffs) ->
    ei_coeff = 0.5 * (a * a for a in euclidean_coeffs).reduce(((b, a) -> b + a), 0)
    @vector [1].concat(euclidean_coeffs).concat([1, ei_coeff])

if typeof module isnt "undefined" and module.exports then module.exports = sph_ga
else window.sph_ga = sph_ga
