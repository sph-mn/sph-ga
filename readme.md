# sph-ga

`sph-ga` is a c library for fundamental calculations in euclidean and conformal geometric algebras.

The implementation is designed around compactness, generality, and directness. It aims to stay structurally simple and adaptable, with explicit data structures and straightforward algebraic operations, rather than pushing complexity into elaborate metaprogramming or other aggressive abstraction machinery.

It provides the core operations of geometric algebra in a compact and explicit form:

- geometric, inner, and exterior products
- euclidean and conformal geometric algebra in the same framework
- sparse multivectors with direct blade-level representation
- reverse, involute, conjugate, inverse, and conformal helpers such as points and rotors

String conversion is included so that multivectors can be written concisely in tests and examples.

Here is a project using the javascript version of `sph-ga` for [n-cube rotation](https://sph.mn/audiovisual/animations/ncubes/main.html).

## build
Build the library and tests:

```sh
./exe/compile
````

## install
Install the library and headers:

```sh
./exe/install
```

Or install as symlinks:

```sh
./exe/install-symlink
```

An optional path prefix can be supplied to install into `$prefix/usr/...` instead of `/usr/...`.

## integrating the library
Include:

```c
#include <sph-ga.h>
```

Then link in the usual way:
```sh
cc your_program.c -lsph-ga
```

## usage
The library is centered around a prepared `ga_t` space instance and explicit `ga_mv_t` multivectors.
## creating a space
A diagonal euclidean metric in three dimensions:
```
ga_t r3;
double metric[3];

metric[0] = 1.0;
metric[1] = 1.0;
metric[2] = 1.0;

r3.dimensions = 0;
ga_init_diagonal((&r3), metric, 3, 0);
```
A conformal space based on the euclidean part `[1, 1, 1]`:
```
ga_t c3;
double metric[3];

metric[0] = 1.0;
metric[1] = 1.0;
metric[2] = 1.0;

c3.dimensions = 0;
ga_init_diagonal((&c3), metric, 3, 1);
```
After use:
```
ga_uninit((&r3));
ga_uninit((&c3));
```
## multivectors
A multivector is represented as:
```
typedef struct {
  ga_blade_id_t id;
  double coeff;
} ga_term_t;

typedef struct {
  ga_term_t *term;
  uint64_t term_count;
} ga_mv_t;
```
Terms are sparse. `id` is a blade bitset.

A caller provides the storage:
```
ga_term_t storage[ga_blade_count];
ga_mv_t mv;

mv.term = storage;
mv.term_count = 0;
```
## simple constructors
Scalar and basis blade macros:
```
ga_term_t storage[ga_blade_count];
ga_mv_t mv;

mv.term = storage;

ga_mv_set_zero((&mv));
ga_mv_set_scalar((&mv), 2.0);
ga_mv_set_basis_blade((&mv), 3, 4.0);
```
Examples:
* scalar `2`
* basis blade `4e3`

## string conversion
String conversion is supported for compact construction, tests, and inspection.
### parse
```
ga_term_t storage[ga_blade_count];
ga_mv_t mv;

mv.term = storage;
ga_mv_from_string((&r3), "1 + e2 + 2e1_3", (&mv));
```
### format
```
char buffer[ga_mv_string_capacity];

ga_mv_to_string((&r3), (&mv), buffer);
```
The string layer is trusted-input / trusted-buffer infrastructure. It is not a defensive parser.
## operations
Binary operations:
* `ga_add`
* `ga_subtract`
* `ga_ep_one`
* `ga_gp_one`
* `ga_ip_one`

Unary operations:
* `ga_reverse`
* `ga_involute`
* `ga_conjugate`

Inverse:
* `ga_inverse`

### example
```
ga_term_t left_storage[ga_blade_count];
ga_term_t right_storage[ga_blade_count];
ga_term_t out_storage[ga_blade_count];
ga_mv_t left;
ga_mv_t right;
ga_mv_t out;

left.term = left_storage;
right.term = right_storage;
out.term = out_storage;

ga_mv_from_string((&r3), "e1", (&left));
ga_mv_from_string((&r3), "e2", (&right));
ga_gp_one((&r3), (&left), (&right), (&out));
```
## conformal helpers
When the space is initialized with `is_conformal = 1`, these helper functions are available:
* `ga_no`
* `ga_ni`
* `ga_point`
* `ga_point_euclidean`
* `ga_rotor`

### `no` and `ni`
```
ga_term_t storage[ga_blade_count];
ga_mv_t mv;
char buffer[ga_mv_string_capacity];

mv.term = storage;

ga_no((&c3), 2.0, (&mv));
ga_mv_to_string((&c3), (&mv), buffer);   /* 2e4 */

ga_ni((&c3), -3.0, (&mv));
ga_mv_to_string((&c3), (&mv), buffer);   /* -3e5 */
```
### conformal point
```
ga_term_t storage[ga_blade_count];
ga_mv_t point;
double euclidean[3];
double extracted[3];
char buffer[ga_mv_string_capacity];

point.term = storage;

euclidean[0] = 2.0;
euclidean[1] = 3.0;
euclidean[2] = 4.0;

ga_point((&c3), euclidean, (&point));
ga_mv_to_string((&c3), (&point), buffer);   /* 2e1 + 3e2 + 4e3 + e4 + 14.5e5 */

ga_point_euclidean((&c3), (&point), extracted);
```
### rotor
For a conformal space based on three euclidean dimensions, the rotor coefficient order is:
* scalar
* `e1_2`
* `e1_3`
* `e2_3`

Example:
```
ga_term_t storage[ga_blade_count];
ga_mv_t rotor;
double coeff[4];
char buffer[ga_mv_string_capacity];

rotor.term = storage;

coeff[0] = 2.0;
coeff[1] = 3.0;
coeff[2] = 4.0;
coeff[3] = 5.0;

ga_rotor((&c3), coeff, (&rotor));
ga_mv_to_string((&c3), (&rotor), buffer);   /* 2 + 3e1_2 + 4e1_3 + 5e2_3 */
```
## status-returning functions
Only functions that can fail return `status_t`.

The current fallible functions are:
* `ga_init_full`
* `ga_init_diagonal`
* `ga_add`
* `ga_subtract`
* `ga_ep_one`
* `ga_gp_one`
* `ga_ip_one`
* `ga_inverse`

Typical pattern:
```
status_t status;
ga_t r3;
double metric[3];

metric[0] = 1.0;
metric[1] = 1.0;
metric[2] = 1.0;

status = ga_init_diagonal((&r3), metric, 3, 0);
if (status.id) {
  return status.id;
}
```
## api
### initialization
```
ga_init_full :: ga metric dimensions is_conformal -> status
ga_init_diagonal :: ga diagonal_metric dimensions is_conformal -> status
ga_uninit :: ga -> void
```
### multivector construction and conversion
```
ga_mv_set_zero :: mv -> macro
ga_mv_set_scalar :: mv value -> macro
ga_mv_set_basis_blade :: mv index value -> macro
ga_mv_from_string :: ga string mv -> void
ga_mv_to_string :: ga mv buffer -> void
```
### algebra
```
ga_add :: ga left right out -> status
ga_subtract :: ga left right out -> status
ga_ep_one :: ga left right out -> status
ga_gp_one :: ga left right out -> status
ga_ip_one :: ga left right out -> status
ga_inverse :: ga in reverse_mv denominator_mv out -> status

ga_reverse :: ga in out -> void
ga_involute :: ga in out -> void
ga_conjugate :: ga in out -> void
```
### conformal helpers
```
ga_no :: ga coeff out -> void
ga_ni :: ga coeff out -> void
ga_point :: ga euclidean_coeff out -> void
ga_point_euclidean :: ga in euclidean_coeff -> void
ga_rotor :: ga coeff out -> void
```
## space properties
`ga_t` exposes these relevant prepared-state fields:
```
dimensions
is_conformal
is_symmetric
is_diagonal
null_vector_count
null_vector_start
id_null
pseudoscalar_id
no_id
ni_id
no_index
ni_index
metric
```
## data model
### basis indices
Basis indices are `1`-based in the human-readable representation:
* `e1`
* `e2`
* `e1_3`

### blade ids
Blade ids are bitsets:
* scalar: `0`
* `e1`: `1`
* `e2`: `2`
* `e1_2`: `3`

### multivectors
Multivectors are sparse lists of `(id, coeff)` terms.
## conformal geometric algebra
With `is_conformal = 1`, the space extends the euclidean part by two null basis vectors appended at the end.

For euclidean dimension `d`, the conformal basis order is:
```
e1, e2, ..., ed, no, ni
```
The default conformal metric added by the library is:
* `ei · ei = 1` for the euclidean part
* `no · no = 0`
* `ni · ni = 0`
* `no · ni = -1`
* `ni · no = -1`

The canonical string form still uses `e...` names. Input also accepts `no` and `ni` as aliases.
## string representation
Examples of valid components:
```
3
20.12
e2
20.12e1
e1_4_6
3e1_4_6
no
ni
```
Examples of valid multivectors:
```
2 + e1 + e2_3
1 + no + -2ni
```
Formatting canonicalizes conformal aliases to `e...` basis names.
## tests
The automated tests are in:
```
src/test/sph-ga.c
```
They include:
* algebra oracle cases for `gp`, `ip`, `ep`
* involution / reverse / conjugate checks
* parser / formatter roundtrip checks
* conformal helper tests

## design
The implementation is intentionally bounded and explicit:
* `ga_max_dimensions` fixes the compilation target
* `ga_accumulator_stack_size` bounds internal algebra scratch depth
* `ga_t` combines declarative space state with runtime scratch state
* concurrent use is achieved by using distinct `ga_t` instances, including copies of prepared instances where desired

## excluded
This library does not provide:
* operator overloading
* code generation
* symbolic simplification
* graphical functions
* defensive parsing of invalid multivector strings
* automatic diagonalization of conformal expressions into another basis

## license
`lgpl3+`
