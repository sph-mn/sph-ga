# sph-ga
this is a javascript library for fundamental calculations in euclidean and conformal geometric algebras.

the focus of this implementation is on flexibility (functional with basic data structures and loose coupling for easy abstraction) as well as compactness and generality (straightforward with no workarounds or limitations).

the collection of automated tests may also be useful for testing other geometric algebra libraries.

**status**
a new approach for calculating the inner product has been implemented and is being tested.

# usage
`compiled/sph_ga.js` contains the javascript version.
use with node.js via `require("./sph_ga.js")` or include the code in html using `<script type="text/javascript" src="sph_ga.js"></script>`.

# usage
the examples are using the [coffeescript](https://coffeescript.org/) javascript syntax.

## creating a space
~~~
r3 = new sph_ga [1, 1, 1]
~~~

this example supplies an euclidean metric. for more options, see the section "customization" below.

## creating vectors
all these types are represented as multivector objects.

`s(number)` creates a scalar.
~~~
s1 = r3.s 1
~~~

`basis(n, scale=1)` creates the basis vector with that index.
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

strings can also be used: `r3.mv_from_string "1 + e2 + 2e1_3"`. this is only for multivector creation and does not support other calculations.

## operations
~~~
e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3

exterior_product = r3.ep e1, e2
geometric_product = r3.gp e1, e2
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

# a multivector of the scalar and one or more 1-blades.
# adds as many basis blades in order as specified.
vector :: (coefficient ...) -> multivector

# create a multivector by specifying multiple blade terms
mv :: (((basis_index ...) coefficient) ...) -> multivector

# calculates the left contraction clifford (geometric) inner product.
# the grade of the first argument must be less than or equal to the grade of the second argument
ip :: multivector ... -> multivector

# calculates the exterior product (also known as the wedge product or outer product)
ep :: multivector ... -> multivector

# calculates the geometric product
gp :: multivector ... -> multivector

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

add :: multivector ... -> multivector
subtract :: multivector ... -> multivector
pseudoscalar :: -> multivector
inverse :: multivector -> multivector
grade :: multivector -> integer
id_from_indices :: (basis_index ...) -> id
id_indices :: (id) -> (basis_index ...)
mv_to_string :: multivector -> string
mv_from_string :: string -> multivector

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
sph_ga has special features for cga. when the conformal option is set to true for a space, it automatically adds two dimensions and uses the conformal metric.

* the null vectors will always be appended to the end of the list of canonical bases: e1, e2, ..., en, no, ni.
* the default conformal metric is a split-signature metric with two off-diagonal -1 terms for the null vector interactions.
  * ei · ei = 1 for 1 <= i <= n
  * no · no = 0
  * ni · ni = 0
  * no · ni = -1
  * ei · n0 = 0 and ei · ni = 0 for 1 <= i <= n
  * ei · ej = 0 for i != j
* the pseudoscalar is defined as: pseudoscalar = e1 * e2 * ... * en * no * ni
  * normalization: pseudoscalar ** 2 = 1   (0 for degenerate metrics, such as with conformal: true)

## usage
~~~
# only the euclidean part of the metric has to be specified.
c3 = new sph_ga [1, 1, 1], conformal: true

# reflecting a point across a plane defined by a normal vector n.
# assume n is a unit vector in the conformal space
e1 = c3.basis 1
e2 = c3.basis 2
n = c3.add e1, e2
point = c3.point 1, 2, 3
# reflection formula: p' = -n * p * n
reflected_point = c3.gp c3.gp(c3.reverse(n), point), n
~~~

## extended api when `conformal` is true
functions
~~~
# create a basis vector for the origin. also known as "n0"
no :: coefficient -> multivector

# create a basis vector for infinity. also known as "n∞"
ni :: coefficient -> multivector

# create a conformal point. requires only the coefficients for the euclidean part
point :: (coefficient ...) -> multivector

# create a rotor. takes the coefficient for the scalar followed by the rotation axes
rotor :: (coefficient ...) -> multivector

# a normalized multivector across all basis blades
normal :: -> multivector
~~~

space properties
~~~
no_index
ni_index
no_id
ni_id
~~~

note that coffeescript does not allow `no` as a variable name.

# customization
## metric tensor
the metric array passed in the options defines the signature of the space. each element represents the square of a basis vector combination:
* `1` for spacelike dimensions
* `-1` for timelike dimensions
* `0` for null dimensions

other values, eg for non-diagonal metrics, are possible. metrics must be symmetric.
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

# string representations
parsing and creating multivector string representations is supported.

valid example components
~~~
20.12
e2
20.12e1
e1_4_6
20.12e1_4_6
~~~

multivector example
~~~
2 + e1 + e2_3
~~~

conversion examples
~~~
mystring = r3.mv_from_string "2 + e1 + e2_3"
r3.mv_to_string mymv
~~~

## ebnf
~~~
input          = [ number ], [ letters ], [ "_" , integer_list ] ;
number         = digit , { digit } , [ "." , digit , { digit } ] ;
letters        = letter , { letter } ;
integer_list   = integer , { "_" , integer } ;
integer        = digit , { digit } ;
digit          = "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" ;
letter         = "a" | "b" | "c" | "d" | ... | "z" ;
~~~

# tests
run via `./exe/tests`. the code for the test cases, data generator, and runner, is located in src/test.coffee.

# excluded
what this library will not provide:
* operator overloading
* code generation
* string notation for complex calculations
* graphical functions

what this library does not do:
* translation of cga expressions into an equivalent minkowski space with a diagonalized metric

# license
lgpl3+

# possible enhancements
* simplify multivector component access
* performance optimizations
* the library could easily be ported to c
