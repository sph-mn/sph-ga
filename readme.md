# sph-ga
this is a javascript library for basic calculations of geometric algebras, including conformal geometric algebra.

the focus of this library is flexibility (functional with basic data structures so that the code can be abstracted easily) and generality (2, 10, or more dimensions should be possible).
two properties which to a certain extend exclude class hierarchies, convenience features, and use-case optimization.

# license
lgpl3+

# usage
compiled/ga.js contains the javascript version.
include the code in the browser using <script > tags or in node.js using require "./ga.js".

# examples
examples using the [coffeescript](https://coffeescript.org/) javascript syntax.

~~~
r3 = new simple_ga [1, 1, 1]
e1 = r3.basis 1
e2 = r3.basis 2
e3 = r3.basis 3
e12 = r3.gp e1, e2
exterior_product_result = r3.ep e1, e2
~~~

~~~
mv_a = [[1, 2], [3, 4]]  # 3e1 + 4e2
mv_b = [[3, 6], [5, -2]]  # 5e1e2 - 2e2e3
inner_product_result = r3.ip mv_a, mv_b
~~~

~~~
[ps_ids, _] = r3.pseudoscalar()
ps_grade = r3.grade ps_ids[0]  # should be 3
ps_squared = r3.gp ps, ps  # should be -1 for this space
~~~

# api
~~~
# the length of the given metric defines the number of dimensions.
constructor :: [integer, ...]:metric -> object

# calculates the geometric product.
gp :: multivector multivector -> multivector

# calculates the inner product as the left contraction.
# the grade of the first argument must be less than or equal to the grade of the second argument.
ip :: multivector multivector -> multivector

# calculates the exterior product. the exterior product is also known as the wedge product or outer product.
ep :: multivector multivector -> multivector

# compute the sandwich product: a * b * a ** -1
sp :: multivector multivector -> multivector

reverse :: multivector -> multivector
involute :: multivector -> multivector
conjugate :: multivector -> multivector
add :: multivector multivector -> multivector
subtract :: multivector multivector -> multivector
basis :: grade -> multivector
blade_id :: grade -> integer
pseudoscalar :: -> multivector
~~~

space object properties
~~~
n: integer  # number of dimensions
metric: [integer, ...]
blade_count: integer
~~~

# data structures
blades are identified by bit-set integers, since this naturally encodes the subset nature of basis vectors. this allows quick access to operations like blade addition, where bitwise or operations can combine different basis vectors or blades. blade identifier 0 is for the scalar.

multivectors are 2-arrays of arrays [blade_identifiers, coefficients]. this is a sparse representation and can be modified to use different number array types like Float64Array.
multivector blade coefficients are therefore accessed via "mv[1][index]".
the compromise of this data structure is that the indices of blade identifiers have to be first found inside the blade_identifier list mv[0] to know the corresponding coefficient index.

# conformal geometric algebra
core geometric algebra operations remain the same in conformal geometric algebra, no special library setup is necessary.

in conformal geometric algebra, the interpretation of core operations and the structure of the space in which they are applied are significantly different due to the conformal models particular choice of metric and dimensional embedding. null vectors square to zero, unlike typical euclidean vectors.

# optimization options
on-demand caching for:
* sign functions for blade grade pairings and given metric
* blade id grades
* blade combinations

other options
* decomposing into lower-dimensional subspaces

# ideas
* the library could easily be ported to c
