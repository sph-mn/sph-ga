class sph_ga
  constructor: (metric) ->
    @n = metric.length
    @metric = metric

  basis: (i) -> if i then [[1 << (i - 1), 1, 1]] else [[0, 1, 0]]
  vector: (coeffs) -> (@basis i for a, i in coeffs)
  blade: (indices, coeff) -> [indices.reduce(((a, i) -> a |= 1 << (i - 1)), 0), coeff, indices.length]
  mv: (terms) -> @blade indices, coeff for [indices, coeff] in terms
  apply_grade_sign: (a, sign_function) -> ([id, coeff * sign_function(grade), grade] for [id, coeff, grade] in a)
  involute: (a) -> @apply_grade_sign a, (grade) -> (-1) ** grade
  conjugate: (a) -> @apply_grade_sign a, (grade) -> (-1) ** ((grade * (grade + 1)) >> 1)
  reverse: (a) -> @apply_grade_sign a, (grade) -> (-1) ** ((grade * (grade - 1)) >> 1)
  add: (a, b) -> @combine a, b, 1
  subtract: (a, b) -> @combine a, b, -1
  sp: (a, b) -> @gp @gp(a, b), @inverse(a)
  pseudoscalar: -> [@blade([1..@n], 1)]

  grade: (a) ->
    n = 0
    while a != 0
      a &= a - 1
      n += 1
    n

  antisymmetric_sign: (a, b) ->
    count = 0
    count_b = 0
    for k in [0...@n - 1]
      count += count_b if (a >> k) & 1
      count_b += 1 if (b >> k) & 1
    if (count & 1) == 0 then 1 else -1

  metric_sign: (a, b) ->
    c = a & b
    sign = 1
    i = 0
    while c != 0
      sign *= @metric[i] if c & 1
      c >>= 1
      i += 1
    sign

  combine: (a, b, scalar = 1) ->
    coeffs = {}
    coeffs[id] = [coeff, grade] for [id, coeff, grade] in a
    for [id, coeff, grade] in b
      if coeffs[id]? then coeffs[id][0] += coeff * scalar
      else coeffs[id] = [coeff * scalar, grade]
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  gp: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        id = id_a ^ id_b
        sign = @antisymmetric_sign(id_a, id_b) * @metric_sign(id_a, id_b)
        coeff = sign * coeff_a * coeff_b
        if coeffs[id]? then coeffs[id][0] += coeff
        else coeffs[id] = [coeff, @grade id]
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  ip: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        if grade_a <= grade_b
          id = id_a ^ id_b
          grade = grade_b - grade_a
          if @grade(id) == grade
            sign = @antisymmetric_sign(id_a, id_b) * @metric_sign(id_a, id_b)
            coeff = sign * coeff_a * coeff_b
            if coeffs[id]? then coeffs[id][0] += coeff
            else coeffs[id] = [coeff, grade]
    [parseInt(id), coeff, grade] for id, [coeff, grade] of coeffs when coeff != 0

  ep: (a, b) ->
    coeffs = {}
    for [id_a, coeff_a, grade_a] in a
      for [id_b, coeff_b, grade_b] in b
        id = id_a ^ id_b
        grade = grade_a + grade_b
        if @grade(id) == grade
          sign = @antisymmetric_sign(id_a, id_b) * @metric_sign(id_a, id_b)
          coeff = sign * coeff_a * coeff_b
          if coeffs[id]? then coeffs[id][0] += coeff
          else coeffs[id] = [coeff, grade]
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
