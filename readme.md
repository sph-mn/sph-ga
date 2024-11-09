# sph-ga
this is a javascript library for fundamental calculations of euclidean and conformal geometric algebra.

the focus of this library is on compactness and flexibility (functional with basic data structures to allow easy abstraction) as well as generality (no conceptual limits on the number of dimensions).

## status
*work in progress*: still in testing.

~~~
Test 1: Reverse: reverse(1) == 1: Success
Test 2: Reverse: reverse(e1) == e1: Success
Test 3: Reverse: reverse(e1e2) == -e1e2: Success
Test 4: Reverse: reverse(e1e2e3) == -e1e2e3: Success
Test 5: Reverse: reverse(1 + e1) == 1 + e1: Success
Test 6: Reverse: reverse(1 + e1 + e1e2) == 1 + e1 - e1e2: Success
Test 7: Reverse: reverse(e1 + e2 + e1e2) == e1 + e2 - e1e2: Success
Test 8: Reverse: reverse(1 + e1 + e2 + e1e2 + e1e3 + e1e2e3) == 1 + e1 + e2 - e1e2 - e1e3 - e1e2e3: Success
Test 9: Reverse Twice: reverse(reverse(e1 + e1e2)) == e1 + e1e2: Success
Test 10: Involute: involute(e1e2e3) == -e1e2e3: Success
Test 11: Conjugate: conjugate(e1e2) == -e1e2: Success
Test 12: ep: Antisymmetry: a ∧ b == - (b ∧ a): Success
Test 13: ep: Antisymmetry: a ∧ a == 0: Success
Test 14: ep: Distributivity: a ∧ (b + c) == a ∧ b + a ∧ c: Success
Test 15: ep: Distributivity: (a + b) ∧ c == a ∧ c + b ∧ c: Success
Test 16: ep: Exterior Product with Scalars: s ∧ A == s A: Success
Test 17: ep: Exterior Product with Scalars: A ∧ s == s A: Success
Test 18: ep: Zero Product with Scalar Zero: 0 ∧ A == 0: Success
Test 19: ep: Zero Product with Scalar Zero: A ∧ 0 == 0: Success
Test 20: ep: Grading: grade(a ∧ Bk) = 1 + k: Success
Test 21: ep: Decomposability: a ∧ b = 0 when a and b are linearly dependent: Success
Test 22: ep: Exterior Product of Higher Grades: e1 ∧ e2 ∧ e3 != 0: Success
Test 23: ep: Exterior Product of Higher Grades: e1 ∧ e1 == 0: Success
Test 24: ep: Exterior Product of Higher Grades: e2 ∧ e2 == 0: Success
Test 25: ep: Exterior Product of Higher Grades: e3 ∧ e3 == 0: Success
Test 26: ep: Exterior product sign conflict test: (e1 ∧ e2) ∧ e3 == -e123: Success
Test 27: ip: Commutativity: a ⋅ b == b ⋅ a: Success
Test 28: ip: Scalar and Multivector: s ⋅ A == 0: Success
Test 29: ip: Multivector and Scalar: A ⋅ s == 0: Success
Test 30: ip: Vector with itself: a ⋅ a == |a|^2: Success
Test 31: ip: Orthogonal Vectors: a ⋅ b == 0: Success
Test 32: ip: Distributivity: a ⋅ (b + c) == a ⋅ b + a ⋅ c: Success
Test 33: ip: Distributivity: (a + b) ⋅ c == a ⋅ c + b ⋅ c: Success
Test 34: ip: Inner Product of Different Grades: A_k ⋅ B_m == 0 when k ≠ m: Success
Test 35: ip: Inner Product of Different Grades: B_m ⋅ A_k == 0 when m ≠ k: Success
Test 36: ip: Basis Vectors Inner Product: e_i ⋅ e_i == 1: Success
Test 37: ip: Basis Vectors Inner Product: e_i ⋅ e_j == 0 when i ≠ j: Success
Test 38: ip: Inner product sign conflict test: Success
Test 39: cga: ep: Anti-Commutativity: e_i ∧ e_j == -e_j ∧ e_i: Success
Test 40: cga: ep: Associativity: (e_i ∧ e_j) ∧ e_k == e_i ∧ (e_j ∧ e_k): Success
Test 41: cga: ep: Distributivity: e_i ∧ (e_j + e_k) == e_i ∧ e_j + e_i ∧ e_k: Success
Test 42: cga: ep: Idempotency: e_i ∧ e_i == 0: Success
Test 43: cga: ep: Exterior product with conformal basis vectors: eo ∧ ei == 1: Success
Test 44: cga: ep: Nested wedge products with conformal vectors creates bivector: eo ∧ e1: Success
Test 45: cga: ep: Nested wedge products with conformal vectors: (eo ∧ e_i) ∧ e_j == eo ∧ (e_i ∧ e_j): Success
Test 46: cga: ep: Higher-Grade Anti-Commutativity: e_i ∧ e_j ∧ e_k == -e_j ∧ e_i ∧ e_k: Success
Test 47: cga: ep: Exterior product of conformal points: P = e0 ∧ e1 ∧ e2 ∧ e3 ∧ eo: Success
Test 48: cga: ip: Inner product of standard basis vectors: e_i · e_j == δ_ij: Success
Test 49: cga: ip: Inner product of same standard basis vectors: e_i · e_i == 1: Success
Test 50: cga: ip: Inner product involving conformal basis vectors: eo · ei == -1: Success
Test 51: cga: ip: Inner product involving conformal basis vectors: eo · eo == 0: Success
Test 52: cga: ip: Inner product involving conformal basis vectors: ei · ei == 0: Success
Test 53: cga: ip: Inner product involving conformal and standard basis vectors: eo · _i == 0: Success
Test 54: cga: ip: Inner product involving conformal and standard basis vectors: ei · e_i == 0: Success
Test 55: cga: ip: Orthogonality of conformal and standard basis vectors: e_i · eo == 0: Success
Test 56: cga: ip: Orthogonality of conformal and standard basis vectors: e_i · ei == 0: Success
Test 57: cga: ip: Inner product of multivectors: e_i · (e_j ∧ e_k) == δ_ij e_k - δ_ik e_j: Success
Test 58: cga: ip: Inner product with eo: eo · (A ∧ ei) == A: Success
Test 59: cga: ip: Inner product with ei: ei · (A ∧ eo) == A: Success
Test 60: cga: ip: Inner product of conformal points: P · eo == 0: Success
Test 61: cga: ip: Inner product of conformal points: P · ei == -0.5 (P · P): Success
Test 62: cga: ip: Scalar inner product: A · B == sum A_i B_i - A+ B- - A- B+: Success
~~~

# license
lgpl3+

# usage
compiled/ga.js contains the javascript version.
use with node.js via `require("./ga.js")` or include the code in html using `<script type="text/javascript" src="ga.js"></script>`.

# examples
examples use the [coffeescript](https://coffeescript.org/) javascript syntax.

## creating basis vectors and using multivector operations
~~~
r3 = new sph_ga [1, 1, 1]

e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3

geometric_product = r3.gp e1, e2
exterior_product = r3.ep e1, e2
~~~

## creating arbitrary multivectors
~~~
# multivector (s e1 e23)
mv1 = r3.mv [
  [[0], 1],  # scalar part: 1
  [[1], 2],  # 2 * e1
  [[2, 3], 3]  # 3 * e23
]

# multivector (e2, e123)
mv2 = r3.mv [
  [[2], 4],  # 4 * e2
  [[1, 2, 3], 5]  # 5 * e123
]

inner_product = r3.ip mv1, mv2
~~~

## accessing components
multivector components are at corresponding array offsets. for instance, mv1 (s e1 e23) has its parts at index 0, 1, and 2 respectively.
the coefficient is always at index 1 of a component.

~~~
scalar_part = mv1[0][1]
e23_coefficient = mv1[2][1]
~~~

## working with the pseudoscalar
~~~
ps = r3.pseudoscalar()
ps_grade = ps[0][2]  # should be 3 for r3
ps_squared = r3.gp ps, ps  # should be -1 for r3
~~~

## accessing space properties
~~~
r3.n  # dimensions
~~~

# api
~~~
# the length of the given metric defines the number of dimensions.
constructor :: array:(integer ...):metric -> object

# create a multivector representing a scalar
s :: number -> multivector

# get a basis vector for index. 1 -> e1, 2 -> e2, ..., n -> en
basis :: integer:index -> multivector

# a multivector based on the bases, for instance s * e1, s * e2, ..., s * en
# it is possible to set only the first n.
vector :: array:(coefficient ...) -> multivector

# create an n-blade. for example, blade([1, 3], 1) creates "e1^e3"
blade :: array:((index ...) number:coeff) -> array:blade:(xor_bits, coefficient, grade)

# create a multivector by specifying multiple blades like for "blade"
mv :: array:((indices coeff) ...) -> multivector

# calculates the geometric product
gp :: multivector multivector -> multivector

# calculates the inner product as the left contraction.
# the grade of the first argument must be less than or equal to the grade of the second argument
ip :: multivector multivector -> multivector

# calculates the exterior product. the exterior product is also known as the wedge product or outer product
ep :: multivector multivector -> multivector

# compute the sandwich product: a * b * a ** -1
sp :: multivector multivector -> multivector

# reverse the order of basis vectors in each blade.
# each blade's coefficient is multiplied by (-1 ** (k * (k - 1) / 2)), where k is the grade.
reverse :: multivector -> multivector

# changes the sign of blades based on their grade.
# each blades coefficient is multiplied by (-1) ** k
involute :: multivector -> multivector

# combines the reverse and involute operations.
# each blades coefficient is multiplied by (-1) ** (k * (k + 1) / 2)
conjugate :: multivector -> multivector

add :: multivector multivector -> multivector
subtract :: multivector multivector -> multivector
pseudoscalar :: -> multivector
inverse :: multivector -> multivector
~~~

space object properties
~~~
n: integer  # number of dimensions
metric: array:(integer ...)
~~~

# notes about the data structures
elements are basic javascript arrays and numbers and there are no restrictions on manipulating them directly.
blades are identified by bit-set integers - one bit for each basis blade combined - since this naturally encodes the subset nature of basis vectors. blade index 0 is for the scalar.
multivectors use a naturally sparse array format. this is however not needed in core operations and it can still be indexed when needed.

the metric array passed to the constructor defines the signature of the space. each element represents the square of a basis vector:
* `1` for spacelike dimensions (basis vector squares to +1)
* `-1` for timelike dimensions (basis vector squares to -1)
* `0` for null dimensions (basis vector squares to 0)

# conformal geometric algebra
in conformal geometric algebra (cga), the interpretation of fundamental operations and the underlying structure of the space are markedly distinct, owing to the conformal model's specific selection of metric and dimensional embedding. although the metric tensor in cga is frequently expressible as a diagonal matrix, such as [1, 1, 1, -1, 0], this representation is inadequate for encapsulating the intricacies inherent to cga due to the incorporation of null vectors. consequently, we define the metric using a complete matrix tensor and explicitly designate the indices corresponding to the null vector components.

## examples
~~~
# initialize the conformal geometric algebra R4,1
c3 = new sph_ga [1, 1, 1, 0, 0], [3, 4]

# create basis vectors e1, e2, e3, e4, and e5
e1 = c3.basis 1
e2 = c3.basis 2
e3 = c3.basis 3
e4 = c3.basis 4  # typically e0
e5 = c3.basis 5  # typically e∞

geometric_product = c3.gp e1, e2
exterior_product = c3.ep e1, e2
inner_product = c3.ip e1, e2

# creating a conformal point
# p = e1 + e2 + e3 + 0.5 * (e4 + e5)
point = c3.add(
  c3.add(e1, e2),
  c3.add(e3, c3.mv [
    [[4], 0.5],  # 0.5 * e4
    [[5], 0.5]   # 0.5 * e5
  ]))

# reflecting a point across a plane defined by a normal vector n.
# assume n is a unit vector in the conformal space
n = c3.add e1, e2
# reflection formula: p' = -n * p * n
reflected_point = c3.gp c3.gp(c3.reverse(n), point), n
~~~

# ideas
* the library could easily be ported to c
* performance optimization
  * caching: the number of blade comparisons and sign computations grows exponentially with each additional dimension. while other libraries often address this by precomputing values tailored to the current space, this approach appears to be effective only up to a finite number of dimensions
  * selective decomposition: decompose computations into lower-dimensional subspaces where possible to reduce complexity
  * conditional short-circuit evaluations: employ conditional checks to bypass unnecessary computations, thereby improving efficiency when certain results can be determined early
