sph_ga = require "./ga.coffee"

array_equal = (a, b) ->
  return false unless Array.isArray(a) and Array.isArray(b)
  return false unless a.length is b.length
  for i in [0...a.length]
    if Array.isArray(a[i]) and Array.isArray(b[i])
      return false unless array_equal(a[i], b[i])
    else
      return false unless a[i] is b[i]
  true

mv_equal = (mv1, mv2) ->
  # Function to compare two multivectors within a specified tolerance
  tolerance = 1e-6
  coeffs1 = {}
  coeffs2 = {}
  for [id, coeff, grade] in mv1
    coeffs1[id.toString()] = coeff
  for [id, coeff, grade] in mv2
    coeffs2[id.toString()] = coeff
  all_ids = {}
  for id of coeffs1
    all_ids[id] = true
  for id of coeffs2
    all_ids[id] = true
  for id of all_ids
    c1 = coeffs1[id] ? 0
    c2 = coeffs2[id] ? 0
    if Math.abs(c1 - c2) > tolerance
      return false
  return true

c3 = new sph_ga [1, 1, 1], conformal: true

run_tests = (tests) ->
  for test, i in tests
    data = if typeof test.data is "function" then test.data() else test.data
    actual =  if typeof test.actual is "function" then test.actual(data) else test.actual
    expected =  if typeof test.expected is "function" then test.expected(data) else test.expected
    is_array = Array.isArray actual
    is_mv = is_array and (!actual.length or (Array.isArray(actual[0]) and 3 == actual[0].length))
    if is_mv then is_equal = mv_equal actual, expected
    else if is_array then is_equal = array_equal actual, expected
    else is_equal = actual == expected
    if is_equal
      console.log "Test #{i + 1}: #{test.title}: Success"
    else
      console.log "Test #{i + 1}: #{test.title}: Failure"
      console.log "  expected", expected if expected
      console.log "  actual  ", c3.mv_string(actual) if actual
      process.exit 1
  console.log "All tests passed."
  process.exit 0

e1 = c3.basis 1
e2 = c3.basis 2
e3 = c3.basis 3

console.log c3

basic_tests = [
  {
    title: "id_bit_indices",
    actual: () -> (c3.id_bit_indices i for i in [0...(2 ** c3.n)])
    expected: [[], [0], [1], [0, 1], [2], [0, 2], [1, 2], [0, 1, 2]]
  }
]

map_grade_factor_tests = [
  {
    title: "Reverse: reverse(1) == 1",
    actual: () -> c3.reverse(c3.vector([1])),
    expected: c3.vector([1])
  }
  {
    title: "Reverse: reverse(e1) == e1",
    actual: () -> c3.reverse(e1),
    expected: e1
  }
  {
    title: "Reverse: reverse(e1e2) == -e1e2",
    actual: () -> c3.reverse(c3.mv([[[1, 2], 1]])),
    expected: c3.mv([[[1, 2], -1]])
  }
  {
    title: "Reverse: reverse(e1e2e3) == -e1e2e3",
    actual: () -> c3.reverse(c3.mv([[[1, 2, 3], 1]])),
    expected: c3.mv([[[1, 2, 3], -1]])
  }
  {
    title: "Reverse: reverse(1 + e1) == 1 + e1",
    actual: () -> c3.reverse(c3.vector([1, 1])),
    expected: c3.vector([1, 1])
  }
  {
    title: "Reverse: reverse(1 + e1 + e1e2) == 1 + e1 - e1e2",
    actual: () -> c3.reverse(c3.mv([[[0], 1], [[1], 1], [[1, 2], 1]])),
    expected: c3.mv([[[0], 1], [[1], 1], [[1, 2], -1]])
  }
  {
    title: "Reverse: reverse(e1 + e2 + e1e2) == e1 + e2 - e1e2",
    actual: () -> c3.reverse(c3.mv([[[1], 1], [[2], 1], [[1, 2], 1]])),
    expected: c3.mv([[[1], 1], [[2], 1], [[1, 2], -1]])
  }
  {
    title: "Reverse: reverse(1 + e1 + e2 + e1e2 + e1e3 + e1e2e3) == 1 + e1 + e2 - e1e2 - e1e3 - e1e2e3",
    actual: () -> c3.reverse(c3.mv([
      [[0], 1],         # Scalar: 1
      [[1], 1],         # e1
      [[2], 1],         # e2
      [[1, 2], 1],      # e1e2
      [[1, 3], 1],      # e1e3
      [[1, 2, 3], 1]    # e1e2e3
    ])),
    expected: c3.mv([
      [[0], 1],         # Scalar: 1
      [[1], 1],         # e1
      [[2], 1],         # e2
      [[1, 2], -1],     # -e1e2
      [[1, 3], -1],     # -e1e3
      [[1, 2, 3], -1]   # -e1e2e3
    ])
  }
  {
    title: "Reverse Twice: reverse(reverse(e1 + e1e2)) == e1 + e1e2",
    actual: () -> c3.reverse(c3.reverse(c3.mv([[[1], 1], [[1, 2], 1]]))),
    expected: c3.mv([[[1], 1], [[1, 2], 1]])
  }
  {
    title: "Involute: involute(e1e2e3) == -e1e2e3",
    actual: () -> c3.involute(c3.mv([[[1, 2, 3], 1]])),
    expected: c3.mv([[[1, 2, 3], -1]])
  }
  {
    title: "Conjugate: conjugate(e1e2) == -e1e2",
    actual: () -> c3.conjugate(c3.mv([[[1, 2], 1]])),
    expected: c3.mv([[[1, 2], -1]])
  }
]

ep_tests = [
  {
    title: "ep: Anticommutativity: a ∧ b == -b ∧ a"
    data: ->
      {
        a: c3.vector([0, 1, 0, 0])  # a = e1
        b: c3.vector([0, 0, 1, 0])  # b = e2
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.negate(c3.ep(d.b, d.a))
  }
  {
    title: "ep: Exterior product of vector with itself is zero: a ∧ a == 0"
    data: ->
      {
        a: c3.vector([0, 1, 2, 3])  # a = e1 + 2e2 + 3e3
      }
    actual: (d) -> c3.ep(d.a, d.a)
    expected: (d) -> c3.s(0)
  }
  {
    title: "ep: Associativity: (a ∧ b) ∧ c == a ∧ (b ∧ c)"
    data: ->
      {
        a: c3.vector([0, 1, 0, 0])  # a = e1
        b: c3.vector([0, 0, 1, 0])  # b = e2
        c: c3.vector([0, 0, 0, 1])  # c = e3
      }
    actual: (d) -> c3.ep(c3.ep(d.a, d.b), d.c)
    expected: (d) -> c3.ep(d.a, c3.ep(d.b, d.c))
  }
  {
    title: "ep: Distributivity over addition: a ∧ (b + c) == a ∧ b + a ∧ c"
    data: ->
      {
        a: c3.vector([0, 1, 0, 0])  # a = e1
        b: c3.vector([0, 0, 1, 0])  # b = e2
        c: c3.vector([0, 0, 0, 1])  # c = e3
      }
    actual: (d) -> c3.ep(d.a, c3.add(d.b, d.c))
    expected: (d) -> c3.add(c3.ep(d.a, d.b), c3.ep(d.a, d.c))
  }
  {
    title: "ep: Scalar multiplication: s ∧ a == s * a"
    data: ->
      {
        s: c3.s(2)
        a: c3.vector([0, 1, 0, 0])  # a = e1
      }
    actual: (d) -> c3.ep(d.s, d.a)
    expected: (d) -> c3.scale(d.a, 2)
  }
  {
    title: "ep: Exterior product of scalars: s ∧ t == 0"
    data: ->
      {
        s: c3.s(2)
        t: c3.s(3)
      }
    actual: (d) -> c3.ep(d.s, d.t)
    expected: (d) -> c3.s(0)
  }
  {
    title: "ep: Exterior product of bivectors: (e1 ∧ e2) ∧ (e2 ∧ e3) == 0"
    data: ->
      {
        b1: c3.mv([[[1,2],1]])  # b1 = e1 ∧ e2
        b2: c3.mv([[[2,3],1]])  # b2 = e2 ∧ e3
      }
    actual: (d) -> c3.ep(d.b1, d.b2)
    expected: (d) -> c3.s(0)
  }
  {
    title: "ep: Exterior product of vector and bivector: a ∧ (e2 ∧ e3)"
    data: ->
      {
        a: c3.vector([0,1,0,0])  # a = e1
        b: c3.mv([[[2,3],1]])  # b = e2 ∧ e3
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.mv([[[1,2,3],1]])  # Expected: e1 ∧ e2 ∧ e3
  }
  {
    title: "ep: Scalar times bivector: s ∧ (e1 ∧ e2) == s * (e1 ∧ e2)"
    data: ->
      {
        s: c3.s(2)
        b: c3.mv([[[1,2],1]])  # b = e1 ∧ e2
      }
    actual: (d) -> c3.ep(d.s, d.b)
    expected: (d) -> c3.scale(d.b, 2)
  }
  {
    title: "ep: Zero result when common basis elements: e1 ∧ e1 == 0"
    data: ->
      {
        e1: c3.basis(1)
      }
    actual: (d) -> c3.ep(d.e1, d.e1)
    expected: (d) -> c3.s(0)
  }
  {
    title: "cga: ep: Exterior product of two scalars is zero"
    data: ->
      {
        s1: c3.s 2
        s2: c3.s 3
      }
    actual: (d) -> c3.ep(d.s1, d.s2)
    expected: (d) -> c3.mv([])
  }
  {
    title: "cga: ep: Exterior product of a vector with itself is zero"
    data: ->
      {
        a: c3.vector [0, 1, 2, 3, 0, 0]  # a = 1e1 + 2e2 + 3e3
      }
    actual: (d) -> c3.ep(d.a, d.a)
    expected: (d) -> c3.mv([])
  }
  {
    title: "cga: ep: Exterior product of e1 and e2 equals e12"
    data: ->
      {
        e1: c3.basis 1
        e2: c3.basis 2
      }
    actual: (d) -> c3.ep(d.e1, d.e2)
    expected: (d) -> c3.mv([[[1, 2], 1]])
  }
  {
    title: "cga: ep: Antisymmetry: a ∧ b == - (b ∧ a)"
    data: ->
      {
        a: c3.basis 1
        b: c3.basis 2
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.scale(c3.ep(d.b, d.a), -1)
  }
  {
    title: "cga: ep: Associativity: a ∧ (b ∧ c) == (a ∧ b) ∧ c"
    data: ->
      {
        a: c3.basis 1
        b: c3.basis 2
        c: c3.basis 3
      }
    actual: (d) -> c3.ep(d.a, c3.ep(d.b, d.c))
    expected: (d) -> c3.ep(c3.ep(d.a, d.b), d.c)
  }
  {
    title: "cga: ep: Exterior product of eo and ei equals eo ∧ ei"
    data: ->
      {
        eo: c3.eo 1
        ei: c3.ei 1
      }
    actual: (d) -> c3.ep(d.eo, d.ei)
    expected: (d) -> c3.mv([[[4, 5], 1]])
  }
  {
    title: "cga: ep: Exterior product of eo with itself is zero"
    data: ->
      {
        eo: c3.basis 4
      }
    actual: (d) -> c3.ep(d.eo, d.eo)
    expected: (d) -> c3.mv([])
  }
  {
    title: "cga: ep: Exterior product of ei with itself is zero"
    data: ->
      {
        ei: c3.basis 5
      }
    actual: (d) -> c3.ep(d.ei, d.ei)
    expected: (d) -> c3.mv([])
  }
  {
    title: "cga: ep: Exterior product of eo and e1 equals eo ∧ e1"
    data: ->
      {
        eo: c3.basis 4
        e1: c3.basis 1
      }
    actual: (d) -> c3.ep(d.eo, d.e1)
    expected: (d) -> c3.mv([[[1, 4], -1]])
  }
  {
    title: "cga: ep: Antisymmetry with null vector: e1 ∧ eo == - (eo ∧ e1)"
    data: ->
      {
        e1: c3.basis 1
        eo: c3.eo(4)
      }
    actual: (d) -> c3.ep(d.e1, d.eo)
    expected: (d) -> c3.mv [[[1, c3.eo_index], 4]]
  }
  {
    title: "cga: ep: Exterior product of a scalar and a vector equals scalar times vector"
    data: ->
      {
        s1: c3.s 2
        v: c3.vector [0, 1, 2, 3, 0, 0]
      }
    actual: (d) -> c3.ep(d.s1, d.v)
    expected: (d) -> c3.scale(d.v, 2)
  }
  {
    title: "cga: ep: Exterior product of a vector and a scalar equals scalar times vector"
    data: ->
      {
        s1: c3.s 2
        v: c3.vector [0, 1, 2, 3]
      }
    actual: (d) -> c3.ep(d.v, d.s1)
    expected: (d) -> c3.scale(d.v, 2)
  }
  {
    title: "cga: ep: Exterior product of a blade with itself is zero"
    data: ->
      {
        a: c3.ep(c3.basis(1), c3.basis(2))
      }
    actual: (d) -> c3.ep(d.a, d.a)
    expected: (d) -> c3.s 0
  }
  {
    title: "cga: ep: Exterior product of linearly dependent vectors is zero"
    data: ->
      {
        a: c3.basis 1
        b: c3.basis(1, 2)
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.s 0
  }
  {
    title: "cga: ep: Exterior product of a vector and a bivector"
    data: ->
      {
        a: c3.basis 1
        b: c3.ep(c3.basis(2), c3.basis(3))
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.mv([[[1, 2, 3], 1]])
  }
  {
    title: "cga: ep: Antisymmetry: a ∧ b == - (b ∧ a) for vector and bivector"
    data: ->
      {
        a: c3.basis 1
        b: c3.ep(c3.basis(3), c3.basis(2))
      }
    actual: (d) -> c3.ep(d.a, d.b)
    expected: (d) -> c3.mv [[[1, 2, 3], -1]]
  }
]

e1 = c3.basis(1, 2)    # 2 * e1
e2 = c3.basis(2, 3)    # 3 * e2
e3 = c3.basis(3, 4)    # 4 * e3
eo = c3.eo(5)          # 5 * eo
ei = c3.ei(6)          # 6 * ei
s2 = c3.s(2)     # scalar 2
s3 = c3.s(3)     # scalar 3
vector_nn = c3.vector([0, 2, 3, 4])  # 2e1 + 3e2 + 4e3
vector_nnn = c3.add(vector_nn, eo)    # vector_nn + 5eo
vector_nnn_ei = c3.add(vector_nn, ei) # vector_nn + 6ei
bivector = c3.ep(e1, e2)              # (2e1) ∧ (3e2)
trivector = c3.ep(bivector, e3)       # bivector ∧ (4e3)

ip_tests = [
  {
    title: "ip: scalar ⋅ scalar"
    actual: -> c3.ip(s2, s3)
    expected: c3.s 0
  }
  {
    title: "ip: scalar ⋅ mv_nn"
    actual: -> c3.ip(s2, e1)
    expected: c3.basis(1, 4)
  }
  {
    title: "ip: mv_nn ⋅ scalar"
    actual: -> c3.ip(e1, s3)
    expected: c3.s 0
  }
  {
    title: "ip: scalar ⋅ mv_n"
    actual: -> c3.ip(s2, eo)
    expected: c3.eo 10
  }
  {
    title: "ip: mv_n ⋅ scalar"
    actual: -> c3.ip(eo, s3)
    expected: c3.s 0
  }
  {
    title: "ip: mv_n ⋅ mv_n"
    actual: -> c3.ip(eo, ei)
    expected: c3.s -30
  }
  {
    title: "ip: mv_n ⋅ mv_nn"
    actual: -> c3.ip(eo, e2)
    expected: c3.s 0
  }
  {
    title: "ip: mv_nn ⋅ mv_n"
    actual: -> c3.ip(e2, eo)
    expected: c3.s 0
  }
  {
    title: "ip: mv_nn ⋅ mv_nn"
    actual: -> c3.ip(e1, e2)
    expected: c3.s 0
  }
  {
    title: "ip: mv_nn ⋅ mv_nnn"
    actual: -> c3.ip(e1, vector_nnn)
    expected: c3.s 4
  }
  {
    title: "ip: mv_nnn ⋅ mv_nn"
    actual: -> c3.ip(vector_nnn, e2)
    expected: c3.s 9
  }
  {
    title: "ip: mv_nnn ⋅ mv_nnn"
    actual: -> c3.ip(vector_nnn, vector_nnn_ei)
    expected: c3.s -1
  }
  {
    title: "ip: mv_n ⋅ mv_nnn"
    actual: -> c3.ip(eo, vector_nn)
    expected: c3.s 0
  }
  {
    title: "ip: mv_nnn ⋅ mv_n"
    actual: -> c3.ip(vector_nn, ei)
    expected: c3.s 0
  }
  {
    title: "ip: bivector ⋅ vector"
    actual: -> c3.ip(bivector, e3)
    expected: c3.s 0
  }
  {
    title: "ip: vector ⋅ bivector"
    actual: -> c3.ip(e3, bivector)
    expected: c3.s 0
  }
  {
    title: "ip: trivector ⋅ bivector"
    actual: -> c3.ip(trivector, bivector)
    expected: c3.s 0
  }
  {
    title: "ip: mv_n ⋅ bivector"
    actual: -> c3.ip(eo, bivector)
    expected: c3.s 0
  }
  {
    title: "ip: bivector ⋅ mv_n"
    actual: -> c3.ip(bivector, ei)
    expected: c3.s 0
  }
]

run_tests [
  ip_tests
  #basic_tests
  #map_grade_factor_tests
  #ep_tests
].flat()
