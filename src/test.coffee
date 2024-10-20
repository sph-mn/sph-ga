sph_ga = require "./ga.coffee"

tolerance = 1e-6

equal = (mv1, mv2) ->
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

test_number = 1

test = (name, condition) ->
  if condition
    console.log "Test #{test_number}: #{name}: Success"
  else
    console.log "Test #{test_number}: #{name}: Failure"
  test_number += 1

# R3 space
r3 = new sph_ga [1, 1, 1]

e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3

# Test e1 * e1 == 1
e1_e1 = r3.gp e1, e1
expected = [[0, 1, 0]]
test "e1 * e1 == 1", equal(e1_e1, expected)

# Test e1 * e2 == e1e2
e1_e2 = r3.gp e1, e2
expected = [[3, 1, 2]]
test "e1 * e2 == e1e2", equal(e1_e2, expected)

# Test e1 * e2 * e1 == -e2
e1_e2_e1 = r3.gp e1_e2, e1
expected = [[2, -1, 1]]
test "e1 * e2 * e1 == -e2", equal(e1_e2_e1, expected)

# Test associativity of geometric product
lhs = r3.gp e1_e2, e3
e2_e3 = r3.gp e2, e3
rhs = r3.gp e1, e2_e3
test "Associativity of geometric product: (e1 * e2) * e3 == e1 * (e2 * e3)", equal(lhs, rhs)

# Test e1 ^ e2 == e1e2
e1_wedge_e2 = r3.ep e1, e2
expected = [[3, 1, 2]]
test "e1 ^ e2 == e1e2", equal(e1_wedge_e2, expected)

# Test e1 | e1 == 1
e1_ip_e1 = r3.ip e1, e1
expected = [[0, 1, 0]]
test "e1 | e1 == 1", equal(e1_ip_e1, expected)

# Test e1 | e2 == 0
e1_ip_e2 = r3.ip e1, e2
expected = []
test "e1 | e2 == 0", equal(e1_ip_e2, expected)

# Test reverse of e1e2 == -e1e2
e1e2_rev = r3.reverse e1_e2
expected = [[3, -1, 2]]
test "Reverse of e1e2 == -e1e2", equal(e1e2_rev, expected)

# Test involute of e1e2e3 == -e1e2e3
e1e2 = r3.gp e1, e2
e1e2e3 = r3.gp e1e2, e3
e1e2e3_inv = r3.involute e1e2e3
expected = [[7, -1, 3]]
test "Involute of e1e2e3 == -e1e2e3", equal(e1e2e3_inv, expected)

# Test conjugate of e1e2 == -e1e2
e1e2_conj = r3.conjugate e1_e2
expected = [[3, -1, 2]]
test "Conjugate of e1e2 == -e1e2", equal(e1e2_conj, expected)

# Test pseudoscalar
ps = r3.pseudoscalar()
expected = [[7, 1, 3]]
test "Pseudoscalar in R3 is e1e2e3", equal(ps, expected)

# Test ps * ps == -1
ps_squared = r3.gp ps, ps
expected = [[0, -1, 0]]
test "ps * ps == -1 in R3", equal(ps_squared, expected)

# Conformal geometric algebra
c3 = new sph_ga [1, 1, 1, 1, -1]

e1 = c3.basis 1
e2 = c3.basis 2
e3 = c3.basis 3
e4 = c3.basis 4
e5 = c3.basis 5

# Create a conformal point: p = e1 + e2 + e3 + 0.5 * (e4 + e5)
point = c3.mv [
  [[1], 1]
  [[2], 1]
  [[3], 1]
  [[4], 0.5]
  [[5], 0.5]
]

# Test point * point == 0
p_squared = c3.gp point, point
expected = []
test "Point squared is zero in conformal space", equal(p_squared, expected)

# Test e4 * e4 == 0
e4_e4 = c3.gp e4, e4
test "e4 * e4 == 0 in conformal space", equal(e4_e4, [])

# Test e5 * e5 == 0
e5_e5 = c3.gp e5, e5
test "e5 * e5 == 0 in conformal space", equal(e5_e5, [])

# Test e4 | e5 == -1
e4_ip_e5 = c3.ip e4, e5
expected = [[0, -1, 0]]
test "e4 | e5 == -1 in conformal space", equal(e4_ip_e5, expected)
