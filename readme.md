# sph-ga

*work in progress*: only works for diagonal metrics so far, not for conformal geometric algebra.

this is a javascript library for basic calculations of geometric algebras, including conformal geometric algebra.

the focus of this library is on compactness and flexibility (functional with basic data structures to allow easy abstraction) as well as generality (no limits on the number of dimensions because of pre-computations). these objectives, to some extent, limit the use of class hierarchies, convenience features, and specific use-case optimizations.

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
core geometric algebra operations remain the same in conformal geometric algebra; no special library setup is necessary.

in conformal geometric algebra, the interpretation of core operations and the structure of the space in which they are applied are significantly different due to the conformal models particular choice of metric and dimensional embedding. null vectors square to zero, unlike typical euclidean vectors.

## examples
~~~
# initialize the conformal geometric algebra R4,1
c3 = new sph_ga [1, 1, 1, 1, -1]

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
* performance optimization through caching. per dimension exponentially increasing counts of blade comparisons and sign calculations are made. other libraries tend to solve this by pre-calculation for the current space, but this seems to only work up to a limited number of dimensions
* performance optimization by selectively decomposing into lower-dimensional subspaces