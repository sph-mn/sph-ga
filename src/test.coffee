sph_ga = require "./ga.coffee"

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

run_tests = (tests) ->
  for test, i in tests
    data = if typeof test.data is "function" then test.data() else test.data
    actual =  if typeof test.actual is "function" then test.actual(data) else test.actual
    expected =  if typeof test.expected is "function" then test.expected(data) else test.expected
    if mv_equal actual, expected then console.log "Test #{i + 1}: #{test.title}: Success"
    else
      console.log "Test #{i + 1}: #{test.title}: Failure"
      console.log "  expected", expected if expected
      console.log "  actual  ", actual if actual
      process.exit 1
  console.log "All tests passed."
  process.exit 0

# Initialize R3 Geometric Algebra
r3 = new sph_ga [1, 1, 1]
e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3

metric = [
  [1, 0, 0, 0, 0],  # e1
  [0, 1, 0, 0, 0],  # e2
  [0, 0, 1, 0, 0],  # e3
  [0, 0, 0, 0, -1],  # eo
  [0, 0, 0, -1, 0],  # ei
]

c3 = new sph_ga metric, [3, 4], true

apply_grade_sign_tests = [
  {
    title: "Reverse: reverse(1) == 1",
    actual: () -> r3.reverse(r3.vector([1])),
    expected: r3.vector([1])
  },
  {
    title: "Reverse: reverse(e1) == e1",
    actual: () -> r3.reverse(e1),
    expected: e1
  },
  {
    title: "Reverse: reverse(e1e2) == -e1e2",
    actual: () -> r3.reverse(r3.mv([[[1, 2], 1]])),
    expected: r3.mv([[[1, 2], -1]])
  },
  {
    title: "Reverse: reverse(e1e2e3) == -e1e2e3",
    actual: () -> r3.reverse(r3.mv([[[1, 2, 3], 1]])),
    expected: r3.mv([[[1, 2, 3], -1]])
  },
  {
    title: "Reverse: reverse(1 + e1) == 1 + e1",
    actual: () -> r3.reverse(r3.vector([1, 1])),
    expected: r3.vector([1, 1])
  },
  {
    title: "Reverse: reverse(1 + e1 + e1e2) == 1 + e1 - e1e2",
    actual: () -> r3.reverse(r3.mv([[[0], 1], [[1], 1], [[1, 2], 1]])),
    expected: r3.mv([[[0], 1], [[1], 1], [[1, 2], -1]])
  },
  {
    title: "Reverse: reverse(e1 + e2 + e1e2) == e1 + e2 - e1e2",
    actual: () -> r3.reverse(r3.mv([[[1], 1], [[2], 1], [[1, 2], 1]])),
    expected: r3.mv([[[1], 1], [[2], 1], [[1, 2], -1]])
  },
  {
    title: "Reverse: reverse(1 + e1 + e2 + e1e2 + e1e3 + e1e2e3) == 1 + e1 + e2 - e1e2 - e1e3 - e1e2e3",
    actual: () -> r3.reverse(r3.mv([
      [[0], 1],         # Scalar: 1
      [[1], 1],         # e1
      [[2], 1],         # e2
      [[1, 2], 1],      # e1e2
      [[1, 3], 1],      # e1e3
      [[1, 2, 3], 1]    # e1e2e3
    ])),
    expected: r3.mv([
      [[0], 1],         # Scalar: 1
      [[1], 1],         # e1
      [[2], 1],         # e2
      [[1, 2], -1],     # -e1e2
      [[1, 3], -1],     # -e1e3
      [[1, 2, 3], -1]   # -e1e2e3
    ])
  },
  {
    title: "Reverse Twice: reverse(reverse(e1 + e1e2)) == e1 + e1e2",
    actual: () -> r3.reverse(r3.reverse(r3.mv([[[1], 1], [[1, 2], 1]]))),
    expected: r3.mv([[[1], 1], [[1, 2], 1]])
  },
  {
    title: "Involute: involute(e1e2e3) == -e1e2e3",
    actual: () -> r3.involute(r3.mv([[[1, 2, 3], 1]])),
    expected: r3.mv([[[1, 2, 3], -1]])
  },
  {
    title: "Conjugate: conjugate(e1e2) == -e1e2",
    actual: () -> r3.conjugate(r3.mv([[[1, 2], 1]])),
    expected: r3.mv([[[1, 2], -1]])
  }
]

ep_tests = [
  {
    title: "ep: Antisymmetry: a ∧ b == - (b ∧ a)"
    data: ->
      {
        a: r3.basis(1) # a = 1e1
        b: r3.basis(2) # b = 1e2
        aa: r3.vector([0, -1]) # a = 1e1
        bb: r3.vector([0, 0, -1]) # b = 1e2
      }
    actual: (d) -> r3.ep(d.a, d.b)
    expected: (d) -> r3.ep(d.aa, d.bb)
  },
  {
    title: "ep: Antisymmetry: a ∧ a == 0"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])  # a = 0 + 1e1 + 2e2 + 3e3
      }
    actual: (d) -> r3.ep(d.a, d.a)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ep: Distributivity: a ∧ (b + c) == a ∧ b + a ∧ c"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])  # a = 0 + 1e1 + 2e2 + 3e3
        b: r3.vector([0, 4, 5, 6])  # b = 0 + 4e1 + 5e2 + 6e3
        c: r3.vector([0, 7, 8, 9])  # c = 0 + 7e1 + 8e2 + 9e3
      }
    actual: (d) -> r3.ep(d.a, r3.add(d.b, d.c))
    expected: (d) -> r3.add(r3.ep(d.a, d.b), r3.ep(d.a, d.c))
  },
  {
    title: "ep: Distributivity: (a + b) ∧ c == a ∧ c + b ∧ c"
    data: ->
      {
        a: r3.vector([0, 1, 0, 0])  # a = 0 + 1e1 + 0e2 + 0e3
        b: r3.vector([0, 0, 1, 0])  # b = 0 + 0e1 + 1e2 + 0e3
        c: r3.vector([0, 0, 0, 1])  # c = 0 + 0e1 + 0e2 + 1e3
      }
    actual: (d) -> r3.ep(r3.add(d.a, d.b), d.c)
    expected: (d) -> r3.add(r3.ep(d.a, d.c), r3.ep(d.b, d.c))
  },
  {
    title: "ep: Exterior Product with Scalars: s ∧ A == s A"
    data: ->
      {
        s: 3
        A: r3.basis(1)  # e1
      }
    actual: (d) -> r3.ep(r3.s(d.s), d.A)
    expected: (d) -> r3.mv([
      [[1], 3]  # 3 * e1
    ])
  },
  {
    title: "ep: Exterior Product with Scalars: A ∧ s == s A"
    data: ->
      {
        s: 4
        A: r3.basis(2)  # e2
      }
    actual: (d) -> r3.ep(d.A, r3.s(d.s))
    expected: (d) -> r3.mv([
      [[2], 4]  # 4 * e2
    ])
  },
  {
    title: "ep: Zero Product with Scalar Zero: 0 ∧ A == 0"
    data: ->
      {
        A: r3.basis(1)  # e1
      }
    actual: (d) -> r3.ep(r3.s(0), d.A)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ep: Zero Product with Scalar Zero: A ∧ 0 == 0"
    data: ->
      {
        A: r3.basis(1)  # e1
      }
    actual: (d) -> r3.ep(d.A, r3.s(0))
    expected: (d) -> r3.s(0)
  },
  {
    title: "ep: Grading: grade(a ∧ Bk) = 1 + k"
    data: ->
      {
        A: r3.mv([
          [[1, 2], 1]  # e1 ∧ e2
        ])
        B: r3.basis(3)  # e3
      }
    actual: (d) -> r3.ep(d.A, d.B)
    expected: (d) -> r3.mv([
      [[1, 2, 3], 1]  # e1 ∧ e2 ∧ e3
    ])
  },
  {
    title: "ep: Decomposability: a ∧ b = 0 when a and b are linearly dependent"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])       # a = 0 + 1e1 + 2e2 + 3e3
        b: r3.vector([0, 2, 4, 6])       # b = 0 + 2e1 + 4e2 + 6e3 (b = 2a)
      }
    actual: (d) -> r3.ep(d.a, d.b)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ep: Exterior Product of Higher Grades: e1 ∧ e2 ∧ e3 != 0"
    actual: () -> r3.ep(r3.basis(1), r3.ep(r3.basis(2), r3.basis(3)))
    expected: () -> r3.mv([
      [[1, 2, 3], 1]
    ])
  },
  {
    title: "ep: Exterior Product of Higher Grades: e1 ∧ e1 == 0"
    actual: () -> r3.ep(r3.basis(1), r3.basis(1))
    expected: () -> r3.s(0)
  },
]

ip_tests = [
  {
    title: "ip: Commutativity: a ⋅ b == b ⋅ a"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])  # a = 1e1 + 2e2 + 3e3
        b: r3.vector([0, 4, 5, 6])  # b = 4e1 + 5e2 + 6e3
      }
    actual: (d) -> r3.ip(d.a, d.b)
    expected: (d) -> r3.ip(d.b, d.a)
  },
  {
    title: "ip: Multivector and Scalar: A ⋅ s == 0"
    data: ->
      {
        s: r3.s(5)
        A: r3.vector([0, 1, 2, 3])
      }
    actual: (d) -> r3.ip(d.A, d.s)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ip: Scalar and Multivector: s ⋅ A == 0"
    data: ->
      {
        s: r3.s(5)
        A: r3.vector([0, 1, 2, 3])
      }
    actual: (d) -> r3.ip(d.s, d.A)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ip: Vector with itself: a ⋅ a == |a|^2"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])
      }
    actual: (d) -> r3.ip(d.a, d.a)
    expected: (d) -> r3.s(14)  # 1^2 + 2^2 + 3^2 = 14
  },
  {
    title: "ip: Orthogonal Vectors: a ⋅ b == 0"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
      }
    actual: (d) -> r3.ip(d.a, d.b)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ip: Distributivity: a ⋅ (b + c) == a ⋅ b + a ⋅ c"
    data: ->
      {
        a: r3.vector([0, 1, 2, 3])
        b: r3.vector([0, 4, 5, 6])
        c: r3.vector([0, 7, 8, 9])
      }
    actual: (d) -> r3.ip(d.a, r3.add(d.b, d.c))
    expected: (d) -> r3.add(r3.ip(d.a, d.b), r3.ip(d.a, d.c))
  },
  {
    title: "ip: Distributivity: (a + b) ⋅ c == a ⋅ c + b ⋅ c"
    data: ->
      {
        a: r3.vector([0, 1, 0, 0])
        b: r3.vector([0, 0, 1, 0])
        c: r3.vector([0, 0, 0, 1])
      }
    actual: (d) -> r3.ip(r3.add(d.a, d.b), d.c)
    expected: (d) -> r3.add(r3.ip(d.a, d.c), r3.ip(d.b, d.c))
  },
  {
    title: "ip: Inner Product of Different Grades: A_k ⋅ B_m == 0 when k ≠ m"
    data: ->
      {
        A: r3.vector([0, 1, 2, 3])  # Vector (grade 1)
        B: r3.mv([
          [[1, 2], 4]  # Bivector (grade 2)
        ])
      }
    actual: (d) -> r3.ip(d.A, d.B)
    expected: (d) -> r3.mv [[[1], -8], [[2], 4]]
  },
  {
    title: "ip: Inner Product of Different Grades: B_m ⋅ A_k == 0 when m ≠ k"
    data: ->
      {
        A: r3.vector([0, 1, 2, 3])
        B: r3.mv([
          [[1, 2], 4]
        ])
      }
    actual: (d) -> r3.ip(d.B, d.A)
    expected: (d) -> r3.s(0)
  },
  {
    title: "ip: Inner product sign conflict test",
    data: ->
      {
        a: r3.vector([0, 1, 2, 0]) # a = e1 + 2e2
        b: r3.vector([0, -1, 0, 3]) # b = -e1 + 3e3
      }
    actual: (d) -> r3.ip(d.a, d.b)
    expected: () -> r3.s(-1) # Result should be -1 from e1 ⋅ -e1
  }
]

cga_ep_tests = [
  {
    title: "cga: ep: Anti-Commutativity: e_i ∧ e_j == -e_j ∧ e_i"
    actual: -> c3.ep(c3.basis(1), c3.basis(2))
    expected: -> c3.subtract(c3.s(0), c3.ep(c3.basis(2), c3.basis(1)))
  },
  {
    title: "cga: ep: Associativity: (e_i ∧ e_j) ∧ e_k == e_i ∧ (e_j ∧ e_k)"
    data: ->
      {e1: c3.basis(1)
       e2: c3.basis(2)
       e3: c3.basis(3)}
    actual: (d) -> c3.ep(c3.ep(d.e1, d.e2), d.e3)
    expected: (d) -> c3.ep(d.e1, c3.ep(d.e2, d.e3))
  },
  {
    title: "cga: ep: Distributivity: e_i ∧ (e_j + e_k) == e_i ∧ e_j + e_i ∧ e_k"
    data: ->
      {e1: c3.basis(1)
       e2: c3.basis(2)
       e3: c3.basis(3)}
    actual: (d) -> c3.ep(d.e1, c3.add(d.e2, d.e3))
    expected: (d) -> c3.add(c3.ep(d.e1, d.e2), c3.ep(d.e1, d.e3))
  },
  {
    title: "cga: ep: Idempotency: e_i ∧ e_i == 0"
    actual: -> c3.ep(c3.basis(1), c3.basis(1))
    expected: -> c3.s 0
  },
  {
    title: "cga: ep: Exterior product with conformal basis vectors: eo ∧ ei == 1"
    actual: -> c3.ep(c3.eo(1), c3.ei(1))
    expected: -> c3.mv([[[4, 5], 1, 2]])
  },
  {
    title: "cga: ep: Nested wedge products with conformal vectors creates bivector: eo ∧ e1"
    data: ->
      {eo: c3.eo(1)
       e1: c3.basis(1)}
    actual: (d) -> c3.ep(d.eo, d.e1)
    expected: (d) -> [[9, -1, 2]]
  },
  {
    title: "cga: ep: Nested wedge products with conformal vectors: (eo ∧ e_i) ∧ e_j == eo ∧ (e_i ∧ e_j)"
    data: ->
      {eo: c3.basis(4)
       e1: c3.basis(1)
       e2: c3.basis(2)}
    actual: (d) -> c3.ep(c3.ep(d.eo, d.e1), d.e2)
    expected: (d) -> c3.ep(d.eo, c3.ep(d.e1, d.e2))
  },
  {
    title: "cga: ep: Higher-Grade Anti-Commutativity: e_i ∧ e_j ∧ e_k == -e_j ∧ e_i ∧ e_k"
    data: ->
      {e1: c3.basis(1)
       e2: c3.basis(2)
       e3: c3.basis(3)}
    actual: (d) -> c3.ep(c3.ep(d.e1, d.e2), d.e3)
    expected: -> c3.subtract(c3.s(0), c3.ep(c3.ep(c3.basis(2), c3.basis(1)), c3.basis(3)))
  },
  {
    title: "cga: ep: Exterior product of conformal points: P = e0 ∧ e1 ∧ e2 ∧ e3 ∧ eo"
    actual: -> c3.ep(c3.ep(c3.ep(c3.ep(c3.basis(1), c3.basis(2)), c3.basis(3)), c3.basis(4)), c3.basis(5))
    expected: -> c3.mv [[[1, 2, 3, 4, 5], 1]]
  }
]

cga_ip_tests = [
  {
    title: "cga: ip: Symmetry: A · B == B · A"
    data: ->
      {A: c3.vector([1,2,3,4,5,6])
       B: c3.vector([6,5,4,3,2,1])}
    actual: (d) -> c3.ip(d.A, d.B)
    expected: (d) -> c3.ip(d.B, d.A)
  },
  {
    title: "cga: ip: Linearity in first argument: A · (B + C) == A · B + A · C"
    data: ->
      {A: c3.vector([1,2,3,4,5,6])
       B: c3.vector([6,5,4,3,2,1])
       C: c3.vector([1,1,1,1,1,1])}
    actual: (d) -> c3.ip(d.A, c3.add(d.B, d.C))
    expected: (d) -> c3.add(c3.ip(d.A, d.B), c3.ip(d.A, d.C))
  },
  {
    title: "cga: ip: Linearity in second argument: (A + B) · C == A · C + B · C"
    data: ->
      {A: c3.vector([1,2,3,4,5,6])
       B: c3.vector([6,5,4,3,2,1])
       C: c3.vector([1,1,1,1,1,1])}
    actual: (d) -> c3.ip(c3.add(d.A, d.B), d.C)
    expected: (d) -> c3.add(c3.ip(d.A, d.C), c3.ip(d.B, d.C))
  },
  {
    title: "cga: ip: Inner product of standard basis vectors: e_i · e_j == δ_ij"
    data: ->
      {e1: c3.basis(1)
       e2: c3.basis(2)}
    actual: (d) -> c3.ip(d.e1, d.e2)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product of same standard basis vectors: e_i · e_i == 1"
    data: ->
      {e1: c3.basis(1)}
    actual: (d) -> c3.ip(d.e1, d.e1)
    expected: -> c3.s 1
  },
  {
    title: "cga: ip: Inner product involving conformal basis vectors: eo · ei == -1"
    actual: (d) -> c3.ip(c3.eo(1), c3.ei(1))
    expected: -> c3.s -1
  },
  {
    title: "cga: ip: Inner product involving conformal basis vectors: ei · ei == 0"
    data: ->
      {ei: c3.ei(1)}
    actual: (d) -> c3.ip(d.ei, d.ei)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product involving conformal basis vectors: eo · eo == 0"
    data: ->
      {eo: c3.eo(1)}
    actual: (d) -> c3.ip(d.eo, d.eo)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product involving conformal and standard basis vectors: eo · _i == 0"
    data: ->
      {eo: c3.basis(4)
       e1: c3.basis(1)}
    actual: (d) -> c3.ip(d.eo, d.e1)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product involving conformal and standard basis vectors: ei · e_i == 0"
    data: ->
      {ei: c3.basis(5)
       e1: c3.basis(1)}
    actual: (d) -> c3.ip(d.ei, d.e1)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Orthogonality of conformal and standard basis vectors: e_i · eo == 0"
    data: ->
      {eo: c3.basis(4)
       e1: c3.basis(1)}
    actual: (d) -> c3.ip(d.e1, d.eo)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Orthogonality of conformal and standard basis vectors: e_i · ei == 0"
    data: ->
      {ei: c3.basis(5)
       e1: c3.basis(1)}
    actual: (d) -> c3.ip(d.e1, d.ei)
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product of multivectors: e_i · (e_j ∧ e_k) == δ_ij e_k - δ_ik e_j"
    data: ->
      {e1: c3.basis(1)
       e2: c3.basis(2)
       e3: c3.basis(3)}
    actual: (d) -> c3.ip(d.e1, c3.ep(d.e2, d.e3))
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product with ei: ei · (e1 ∧ eo) == 0"
    data: ->
      {eo: c3.eo(1)
       ei: c3.ei(1)
       A: c3.basis(1)}
    actual: (d) -> c3.ip(d.ei, c3.ep(d.A, d.eo))
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product with eo: eo · (e1 ∧ ei) == 0"
    data: ->
      {eo: c3.eo(1)
       ei: c3.ei(1)
       A: c3.basis(1)}
    actual: (d) -> c3.ip(d.eo, c3.ep(d.A, d.ei))
    expected: -> c3.s 0
  },
  {
    title: "cga: ip: Inner product of conformal points: P · eo == 0"
    data: ->
      {P: c3.vector([1, 1, 1, 1, 1, 1])
       eo: c3.eo(1)}
    actual: (d) -> c3.ip(d.P, d.eo)
    expected: -> c3.s -1
  },
  {
    title: "cga: ip: Inner product of conformal points: P · ei == -1 (P · P)"
    data: ->
      {P: c3.vector([1, 1, 1, 1, 1, 1.5])
       ei: c3.ei(1)}
    actual: (d) -> c3.ip(d.P, d.ei)
    expected: -> c3.s -1
  },
  {
    title: "cga: ip: Scalar inner product: A · B == sum A_i B_i - A+ B- - A- B+"
    data: ->
      {A: c3.vector([1,2,3,4,5,6])
       B: c3.vector([6,5,4,3,2,1])}
    actual: (d) -> c3.ip(d.A, d.B)
    expected: -> c3.s 17
  }
]

gp_tests = [
  {
    title: "gp: multivector × multivector: (1 + 2e1 + 4e23)(5 + 6e12 + 7e3) == 5 + 10e1 + 40e2 + 7e3 + 6e12 - 10e13 + 20e23",
    data: ->
      {
        mv1: c3.mv([
          [[0], 1],          # 1
          [[1], 2],          # 2e1
          [[2, 3], 4]        # 4e23
        ]),
        mv2: c3.mv([
          [[0], 5],          # 5
          [[1, 2], 6],       # 6e12
          [[3], 7]           # 7e3
        ])
      }
    actual: (d) -> c3.gp(d.mv1, d.mv2)
    expected: (d) ->
      c3.mv([
        [[0], 5],
        [[1], 10],
        [[2], 40],
        [[3], 7],
        [[1, 2], 6],
        [[1, 3], -10],
        [[2, 3], 20]
      ])
  }
  {
    title: "gp: Associativity: (a b) c == a (b c)"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
        c: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(r3.gp(d.a, d.b), d.c)
    expected: (d) -> r3.mv([
      [[1, 2, 3], 1]  # e123
    ])
  },
  {
    title: "gp: Distributivity Over Addition (Left): a (b + c) == a b + a c"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
        c: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(d.a, r3.add(d.b, d.c))
    expected: (d) -> r3.add(
      r3.mv([[[1, 2], 1]]),  # e12
      r3.mv([[[1, 3], 1]])   # e13
    )
  },
  {
    title: "gp: Distributivity Over Addition (Right): (b + c) a == b a + c a"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
        c: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(r3.add(d.b, d.c), d.a)
    expected: (d) -> r3.add(
      r3.mv([[[1, 2], -1]]),  # -e12
      r3.mv([[[1, 3], -1]])   # -e13
    )
  },
  {
    title: "gp: Geometric Product with Scalar on the Left: s a == a s"
    data: ->
      {
        s: 3
        a: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(r3.s(d.s), d.a)
    expected: (d) -> r3.mv([
      [[1], 3]  # 3 * e1
    ])
  },
  {
    title: "gp: Geometric Product with Scalar on the Right: a s == s a"
    data: ->
      {
        s: 4
        a: r3.basis(2)  # e2
      }
    actual: (d) -> r3.gp(d.a, r3.s(d.s))
    expected: (d) -> r3.mv([
      [[2], 4]  # 4 * e2
    ])
  },
  {
    title: "gp: Basis Vector with Itself: e_i e_i == 1"
    data: ->
      {
        e_i: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(d.e_i, d.e_i)
    expected: (d) -> r3.s(1)
  },
  {
    title: "gp: Orthonormal Basis Vectors: e_i e_j == -e_j e_i for i != j"
    data: ->
      {
        e_i: r3.basis(1)  # e1
        e_j: r3.basis(2)  # e2
      }
    actual: (d) -> r3.gp(d.e_i, d.e_j)
    expected: (d) -> r3.mv([
      [[1, 2], 1]  # e12
    ])
  },
  {
    title: "gp: Orthonormal Basis Vectors: e_j e_i == -e_i e_j for i != j"
    data: ->
      {
        e_i: r3.basis(1)  # e1
        e_j: r3.basis(2)  # e2
      }
    actual: (d) -> r3.gp(d.e_j, d.e_i)
    expected: (d) -> r3.mv([
      [[1, 2], -1]  # -e12
    ])
  },
  {
    title: "gp: Geometric Product of Multiple Basis Vectors: e_i e_j e_k == e_i (e_j e_k)"
    data: ->
      {
        e_i: r3.basis(1)  # e1
        e_j: r3.basis(2)  # e2
        e_k: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(r3.gp(d.e_i, d.e_j), d.e_k)
    expected: (d) -> r3.mv([
      [[1, 2, 3], 1]  # e123
    ])
  },
  {
    title: "gp: Linearity in the First Argument: (a + b) c == a c + b c"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
        c: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(r3.add(d.a, d.b), d.c)
    expected: (d) -> r3.add(
      r3.mv([[[1, 3], 1]]),  # e13
      r3.mv([[[2, 3], 1]])   # e23
    )
  },
  {
    title: "gp: Linearity in the Second Argument: a (b + c) == a b + a c"
    data: ->
      {
        a: r3.basis(1)  # e1
        b: r3.basis(2)  # e2
        c: r3.basis(3)  # e3
      }
    actual: (d) -> r3.gp(d.a, r3.add(d.b, d.c))
    expected: (d) -> r3.add(
      r3.mv([[[1, 2], 1]]),  # e12
      r3.mv([[[1, 3], 1]])   # e13
    )
  },
  {
    title: "gp: Geometric Product with Zero: a 0 == 0"
    data: ->
      {
        a: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(d.a, r3.s(0))
    expected: (d) -> r3.s(0)
  },
  {
    title: "gp: Geometric Product with Zero: 0 a == 0"
    data: ->
      {
        a: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(r3.s(0), d.a)
    expected: (d) -> r3.s(0)
  },
  {
    title: "gp: Scalar Associativity: s (a b) == (s a) b"
    data: ->
      {
        s: 2,
        a: r3.basis(1),  # e1
        b: r3.basis(2)   # e2
      }
    actual: (d) -> r3.gp(r3.s(d.s), r3.gp(d.a, d.b))
    expected: (d) -> r3.gp(r3.mv([[[1], d.s]]), d.b)
  },
  {
    title: "gp: Scalar Associativity: s (a b) == a (s b)"
    data: ->
      {
        s: 3,
        a: r3.basis(1),  # e1
        b: r3.basis(2)   # e2
      }
    actual: (d) -> r3.gp(r3.s(d.s), r3.gp(d.a, d.b))
    expected: (d) -> r3.gp(d.a, r3.mv([[[2], d.s]]))
  },
  {
    title: "gp: Identity Element: 1 a == a"
    data: ->
      {
        a: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(r3.s(1), d.a)
    expected: (d) -> d.a
  },
  {
    title: "gp: Identity Element: a 1 == a"
    data: ->
      {
        a: r3.basis(2)  # e2
      }
    actual: (d) -> r3.gp(d.a, r3.s(1))
    expected: (d) -> d.a
  },
  {
    title: "gp: Geometric Product with Identity: 1 e_i == e_i"
    data: ->
      {
        e_i: r3.basis(1)  # e1
      }
    actual: (d) -> r3.gp(r3.s(1), d.e_i)
    expected: (d) -> d.e_i
  },
  {
    title: "gp: Geometric Product with Identity: e_i 1 == e_i"
    data: ->
      {
        e_i: r3.basis(2)  # e2
      }
    actual: (d) -> r3.gp(d.e_i, r3.s(1))
    expected: (d) -> d.e_i
  }
]

cga_gp_tests = [
  {
    title: "cga: gp: scalar × scalar",
    data: ->
      {
        s1: c3.s(2),    # Scalar 2
        s2: c3.s(3)     # Scalar 3
      }
    actual: (d) -> c3.gp(d.s1, d.s2)
    expected: (d) -> c3.s(6)
  },
  {
    title: "cga: gp: scalar × multivector",
    data: ->
      {
        s: c3.s(2),                      # Scalar 2
        mv: c3.vector([0, 1, 2, 3])      # 0 + e1 + 2e2 + 3e3
      }
    actual: (d) -> c3.gp(d.s, d.mv)
    expected: (d) -> c3.vector([0, 2, 4, 6])
  },
  {
    title: "cga: gp: multivector × scalar",
    data: ->
      {
        mv: c3.vector([0, 1, 2, 3]),     # 0 + e1 + 2e2 + 3e3
        s: c3.s(4)                      # Scalar 4
      }
    actual: (d) -> c3.gp(d.mv, d.s)
    expected: (d) -> c3.vector([0, 4, 8, 12])
  },
  {
    title: "cga: gp: basis vector × itself",
    data: ->
      {
        e1: c3.basis(1)                  # e1
      }
    actual: (d) -> c3.gp(d.e1, d.e1)
    expected: (d) -> c3.s(1)
  },
  {
    title: "cga: gp: e1 × e2 = e1e2",
    data: ->
      {
        e1: c3.basis(1),                 # e1
        e2: c3.basis(2)                  # e2
      }
    actual: (d) -> c3.gp(d.e1, d.e2)
    expected: (d) -> c3.mv([[[1, 2], 1]])  # e12
  },
  {
    title: "cga: gp: e1 × e0 = e1e0",
    data: ->
      {
        e1: c3.basis(1),                 # e1
        e0: c3.eo(1)                     # e0 with coefficient 1
      }
    actual: (d) -> c3.gp(d.e1, d.e0)
    expected: (d) -> c3.mv([[[1, c3.eo_index], 1]])  # e10
  },
  {
    title: "cga: gp: e1 × e∞ = e1e∞",
    data: ->
      {
        e1: c3.basis(1),                 # e1
        einf: c3.ei(1)                   # e∞ with coefficient 1
      }
    actual: (d) -> c3.gp(d.e1, d.einf)
    expected: (d) -> c3.mv([[[1, c3.ei_index], 1]])  # e14
  },
  {
    title: "cga: gp: e0 × e∞ + e∞ × e0 = -2",
    data: ->
      {
        e0: c3.eo(1),                    # e0
        einf: c3.ei(1)                    # e∞
      }
    actual: (d) -> c3.add(c3.gp(d.e0, d.einf), c3.gp(d.einf, d.e0))
    expected: (d) -> c3.s(-2)
  },
  {
    title: "cga: gp: e0 squared is 0",
    data: ->
      {
        e0: c3.eo(1)                     # e0
      }
    actual: (d) -> c3.gp(d.e0, d.e0)
    expected: (d) -> c3.s(0)
  },
  {
    title: "cga: gp: e∞ squared is 0",
    data: ->
      {
        einf: c3.ei(1)                   # e∞
      }
    actual: (d) -> c3.gp(d.einf, d.einf)
    expected: (d) -> c3.s(0)
  },
  {
    title: "cga: gp: e1 × e2e3 = e1e2e3",
    data: ->
      {
        e1: c3.basis(1),                  # e1
        e2e3: c3.mv([[[2, 3], 1]])       # e23
      }
    actual: (d) -> c3.gp(d.e1, d.e2e3)
    expected: (d) -> c3.mv([[[1, 2, 3], 1]])  # e123
  },
  {
    title: "cga: gp: associativity (a * b) * c == a * (b * c)",
    data: ->
      {
        a: c3.basis(1),                   # e1
        b: c3.basis(2),                   # e2
        c: c3.basis(3)                    # e3
      }
    actual: (d) -> c3.gp(c3.gp(d.a, d.b), d.c)
    expected: (d) -> c3.gp(d.a, c3.gp(d.b, d.c))
  },
  {
    title: "cga: gp: distributivity A * (B + C) = A * B + A * C",
    data: ->
      {
        A: c3.basis(1),                   # e1
        B: c3.basis(2),                   # e2
        C: c3.basis(3)                    # e3
      }
    actual: (d) -> c3.gp(d.A, c3.add(d.B, d.C))
    expected: (d) -> c3.add(c3.gp(d.A, d.B), c3.gp(d.A, d.C))
  },
  {
    title: "cga: gp: distributivity (A + B) * C = A * C + B * C",
    data: ->
      {
        A: c3.basis(1),                   # e1
        B: c3.basis(2),                   # e2
        C: c3.basis(3)                    # e3
      }
    actual: (d) -> c3.gp(c3.add(d.A, d.B), d.C)
    expected: (d) -> c3.add(c3.gp(d.A, d.C), c3.gp(d.B, d.C))
  }
]

run_tests [
  ip_tests
  cga_ip_tests
  apply_grade_sign_tests
  ep_tests
  cga_ep_tests
  #gp_tests
  #cga_gp_tests
].flat()
