# sph-ga
this is a javascript library for fundamental calculations of euclidean and conformal geometric algebras. custom algebras are possible.

the focus of this implementation is on compactness and flexibility (functional, with basic data structures to allow easy abstraction) as well as compactness and generality (straightforward and no limitations).

the collection of automated test might also be useful for testing other geometric algebra libraries.

## status
*work in progress*: still in development.

# license
lgpl3+

# usage
compiled/ga.js contains the javascript version.
use with node.js via `require("./ga.js")` or include the code in html using `<script type="text/javascript" src="ga.js"></script>`.

# usage
the examples are using the [coffeescript](https://coffeescript.org/) javascript syntax.

## creating a space
~~~
r3 = new sph_ga [1, 1, 1]
~~~

this example supplies an euclidean (diagonal) metric. for more options, see the section "customization" below.

## creating vectors
all types are represented as multivector objects.

`s(number)` creates a scalar.
~~~
s1 = r3.s 1
~~~

`basis(n, coefficient=1)` creates the basis vector with that index. index 0 is the scalar.
~~~
e3 = r3.basis 3
~~~

`vector([coefficient, ...])` creates a multivector for as many bases in order as specified. the first array element is the scalar.
~~~
a = r3.vector [0, 2, 3, 4]   # a = 0 + 2e1 + 3e2 + 4e3
~~~

`mv([[[basis_index ...], coefficient], ...])` creates custom multivectors where multiple exterior-product-combined basis indices can be specified for each blade. blades can be specified arbitrarily in any order.
~~~
a = r3.mv [
  [[0], 1]  # scalar part: 1
  [[2], 3]  # 3 * e2
  [[2, 3], 4]  # 4 * e23
]
~~~

## operations
~~~
e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3

geometric_product = r3.gp e1, e2
exterior_product = r3.ep e1, e2
~~~

~~~
# 1 + 2 * e1 + 3 * e23
a = r3.mv [
  [[0], 1]
  [[1], 2]
  [[2, 3], 3]
]

# 4 * e2 + 5 * e123
b = r3.mv [
  [[2], 4]
  [[1, 2, 3], 5]
]

inner_product = r3.ip a, b
~~~

see section "api" for all available functions.

## accessing components
blades are accesible by `id` in a vector object.

`get(multivector, id)` extracts a blade with a specific id or returns `null` if it is not included.

~~~
a = r3.s -4
scalar_coefficient = blade_coeff get a, 0  # -4
~~~

## working with the pseudoscalar
~~~
ps = r3.pseudoscalar()
ps_grade = r3.grade ps  # should be 3 for r3
ps_squared = r3.gp ps, ps  # should be -1 for r3
~~~

## accessing space properties
~~~
r3.n  # dimensions
~~~

## data types
sph_ga uses three compound types: space (object), multivector (array[][3]), and blade (array[3]).

# api
in this library, round bracket enclosed lists are javascript arrays. note that the functions do not validate input arguments and expect strict adherence to the [contract](https://en.wikipedia.org/wiki/Design_by_contract). the listing uses this [type signature format](https://sph.mn/computer/designs/type-signature.html).

~~~
# the number of dimensions is defined by the length of the metric.
constructor :: metric options -> object

# create a scalar
s :: coefficient -> multivector

# get a basis vector of index. 1 -> e1, 2 -> e2, ..., n -> en
basis :: basis_index -> multivector

# a multivector of 1-blades, for instance s * e1, s * e2, ..., s * en
# adds as many basis blades in order as specified.
vector :: (coefficient ...) -> multivector

# create a multivector by specifying multiple blades like for
mv :: (((basis_index ...) coefficient) ...) -> multivector

# calculates the inner product using the left contraction rule.
# the grade of the first argument must be less than or equal to the grade of the second argument
ip :: multivector multivector -> multivector

# calculates the exterior product (also known as the wedge product or outer product)
ep :: multivector multivector -> multivector

# calculates the geometric product
gp :: multivector multivector -> multivector

# compute the sandwich product: a * b * a ** -1
sp :: multivector multivector -> multivector

# reverse the order of basis vectors of each blade.
# each coefficient is multiplied by (-1 ** (k * (k - 1) / 2)), where k is the grade
reverse :: multivector -> multivector

# changes the sign of blades based on their grade.
# each coefficient is multiplied by (-1) ** k
involute :: multivector -> multivector

# combines the reverse and involute operations.
# each coefficient is multiplied by (-1) ** (k * (k + 1) / 2)
conjugate :: multivector -> multivector

add :: multivector multivector -> multivector
subtract :: multivector multivector -> multivector
pseudoscalar :: -> multivector
inverse :: multivector -> multivector
grade :: multivector -> integer
id_from_indices :: (basis_index ...) -> id

# accessors
get :: multivector id -> blade/null
blade_id :: blade -> id
blade_coeff :: blade -> coefficient
blade_grade :: blade -> integer
~~~

space object properties
~~~
n: integer  # number of dimensions
metric: ((integer) ...)  # always an n * n array
pseudoscalar_id
~~~

# data types and structures
~~~
basis_index: integer
id: integer:bitset
coefficient: number
blade: array:(id coefficient grade)
multivector: array:(blade ...)
~~~

# conformal geometric algebra
sph_ga has special features for cga.

## usage
~~~
# initialize the conformal geometric algebra R4,1.
# only the euclidean part of the metric has to be specified.
c3 = new sph_ga [1, 1, 1], conformal: true

# create basis vectors e1, e2, e3, e4, and e5
e1 = c3.basis 1
e2 = c3.basis 2
e3 = c3.basis 3
eo = c3.eo 1
ei = c3.ei 1

# creating a conformal point. p = e1 + e2 + e3 + 0.5 * (e4 + e5)
point = c3.add(
  c3.add(e1, e2),
  c3.add(e3, c3.mv [
    [[c3.eo_index], 0.5],  # 0.5 * e4
    [[c3.ei_index], 0.5]   # 0.5 * e5]))

# reflecting a point across a plane defined by a normal vector n.
# assume n is a unit vector in the conformal space
n = c3.add e1, e2
# reflection formula: p' = -n * p * n
reflected_point = c3.gp c3.gp(c3.reverse(n), point), n
~~~

## extended api when `conformal` is true
functions
~~~
# create a basis vector for the origin. also known as "e0" and "e+"
eo :: coefficient -> multivector

# create a basis vector for infinity. also known as "e∞"" and "e-"
ei :: coefficient -> multivector
~~~

space properties
~~~
eo_index
ei_index
eo_id
ei_id
~~~

# customization

## metric tensor
the metric array passed in the options defines the signature of the space. each element represents the square of a basis vector combination:
* `1` for spacelike dimensions
* `-1` for timelike dimensions
* `0` for null dimensions

other values, eg for non-orthogonal metrics, are possible.
diagonal metrics can be configured using flat arrays. for example, [1, 1, 1].
custom tensors can be provided as an "n * n" array. example for five dimensions:
~~~
metric = [
  [1, 0, 0, 0, 0],
  [0, 1, 0, 0, 0],
  [0, 0, 1, 0, 0],
  [0, 0, 0, 0, -1],
  [0, 0, 0, -1, 0],
]

space = new sph_ga metric
~~~

# ideas
* the library could easily be ported to c
* more performance optimizations. possibly decompose computations into lower-dimensional subspaces
