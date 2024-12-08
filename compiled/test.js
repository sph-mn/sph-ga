// Generated by CoffeeScript 2.7.0
var basic_tests, c3, ep_tests, gen, ip_tests, map_grade_factor_tests, runner, sph_ga, test_data_generator_class, test_runner_class;

sph_ga = require("./ga.coffee");

test_runner_class = class test_runner_class {
  constructor(space) {
    this.space = space;
  }

  array_equal(a, b) {
    var i, j, ref;
    if (!(Array.isArray(a) && Array.isArray(b))) {
      return false;
    }
    if (a.length !== b.length) {
      return false;
    }
    for (i = j = 0, ref = a.length; (0 <= ref ? j < ref : j > ref); i = 0 <= ref ? ++j : --j) {
      if (Array.isArray(a[i]) && Array.isArray(b[i])) {
        if (!array_equal(a[i], b[i])) {
          return false;
        }
      } else {
        if (a[i] !== b[i]) {
          return false;
        }
      }
    }
    return true;
  }

  mv_equal(mv1, mv2) {
    var all_ids, c1, c2, coeff, coeffs1, coeffs2, grade, id, j, l, len, len1, ref, ref1, tolerance;
    // Function to compare two multivectors within a specified tolerance
    tolerance = 1e-6;
    coeffs1 = {};
    coeffs2 = {};
    for (j = 0, len = mv1.length; j < len; j++) {
      [id, coeff, grade] = mv1[j];
      coeffs1[id.toString()] = coeff;
    }
    for (l = 0, len1 = mv2.length; l < len1; l++) {
      [id, coeff, grade] = mv2[l];
      coeffs2[id.toString()] = coeff;
    }
    all_ids = {};
    for (id in coeffs1) {
      all_ids[id] = true;
    }
    for (id in coeffs2) {
      all_ids[id] = true;
    }
    for (id in all_ids) {
      c1 = (ref = coeffs1[id]) != null ? ref : 0;
      c2 = (ref1 = coeffs2[id]) != null ? ref1 : 0;
      if (Math.abs(c1 - c2) > tolerance) {
        return false;
      }
    }
    return true;
  }

  run(...tests) {
    var actual, data, expected, i, is_array, is_equal, is_mv, j, len, ref, test;
    ref = tests.flat();
    for (i = j = 0, len = ref.length; j < len; i = ++j) {
      test = ref[i];
      data = typeof test.data === "function" ? test.data() : test.data;
      actual = typeof test.actual === "function" ? test.actual(data) : test.actual;
      expected = typeof test.expected === "function" ? test.expected(data) : test.expected;
      is_array = Array.isArray(actual);
      is_mv = is_array && (!actual.length || (Array.isArray(actual[0]) && 3 === actual[0].length));
      if (is_mv) {
        is_equal = mv_equal(actual, expected);
      } else if (is_array) {
        is_equal = array_equal(actual, expected);
      } else {
        is_equal = actual === expected;
      }
      if (is_equal) {
        console.log(`Test ${i + 1}: ${test.title}: Success`);
      } else {
        console.log(`Test ${i + 1}: ${test.title}: Failure`);
        if (expected) {
          console.log("  expected", expected);
        }
        if (actual) {
          console.log("  actual  ", this.space.mv_string(actual));
        }
        process.exit(1);
      }
    }
    console.log("All tests passed.");
    return process.exit(0);
  }

};

test_data_generator_class = class test_data_generator_class {
  random_integer(min, max) {
    return Math.floor(Math.random() * (max - min + 1)) + min;
  }

  random_element(a) {
    return a[Math.floor(Math.random() * a.length / 1)];
  }

  shuffle(array) {
    return array.sort(function() {
      return Math.random() - 0.5;
    });
  }

  id_type(id) {
    var type;
    if (id & ~this.space.id_null) {
      if (id & this.space.id_null) {
        return type = "en";
      } else {
        return type = "e";
      }
    } else {
      if (id & this.space.id_null) {
        return type = "n";
      } else {
        return type = "s";
      }
    }
  }

  constructor(space) {
    var a, b, g, id, j, k, ref, t, v;
    this.cycled_coeff = this.cycled_sequence([2, 3, 4, 5, 6, 7, 8, 9]);
    this.space = space;
    this.n = space.n;
    this.blade_ids = (function() {
      var j, ref, results;
      results = [];
      for (j = 0, ref = this.n; (0 <= ref ? j <= ref : j >= ref); 0 <= ref ? j++ : j--) {
        results.push({
          s: [],
          e: [],
          n: [],
          en: []
        });
      }
      return results;
    }).call(this);
    for (id = j = 0, ref = 2 ** this.n; (0 <= ref ? j < ref : j > ref); id = 0 <= ref ? ++j : --j) {
      g = space.id_grade(id);
      t = this.id_type(id);
      this.blade_ids[g][t].push(id);
    }
    this.blade_id_next = (function() {
      var l, len, ref1, results;
      ref1 = this.blade_ids;
      results = [];
      for (l = 0, len = ref1.length; l < len; l++) {
        a = ref1[l];
        b = {};
        for (k in a) {
          v = a[k];
          b[k] = this.cycled_sequence(v);
        }
        results.push(b);
      }
      return results;
    }).call(this);
  }

  cycled_sequence(a) {
    var i;
    i = 0;
    return function() {
      var b;
      b = a[i];
      i += 1;
      i = i < a.length ? i : 0;
      return b;
    };
  }

  blade(grade, type, overlap_id) {
    var id, ids, ids_with_overlap;
    // type
    //   s: scalar
    //   e: euclidean
    //   n: null
    //   en: mixed
    if (overlap_id) {
      ids = this.blade_ids[grade][type];
      ids_with_overlap = ids.filter(function(a) {
        return a & overlap_id;
      });
      id = this.random_element(ids_with_overlap);
    } else {
      id = this.blade_id_next[grade][type]();
    }
    if (id == null) {
      return null;
    }
    return [this.space.id_indices(id), this.cycled_coeff()];
  }

  derived_blade(blade_a, grade_b, type_b, overlap_type) {
    var blade_b, grade_a, id_a, type_a;
    // overlap_type
    //   0: unspecified overlap
    //   1: partial overlap
    //   2: full overlap
    id_a = this.space.id_from_indices(blade_a[0]);
    type_a = this.id_type(id_a);
    grade_a = this.space.id_grade(id_a);
    switch (overlap_type) {
      case 0:
        return this.blade(grade_b, type_b);
      case 1:
        return this.blade(grade_b, type_b, id_a);
      case 2:
        if (!(grade_a === grade_b && type_a === type_b)) {
          return null;
        }
        blade_b = Array.from(blade_a);
        blade_b[1] = this.cycled_coeff();
        return blade_b;
    }
  }

  mv(config) {
    var a, first, previous, rest;
    // [[grade, type, overlap_type], ...]
    first = this.blade.apply(this, config[0].slice(0, 2));
    rest = config.slice(1);
    previous = first;
    rest = (function() {
      var j, len, results;
      results = [];
      for (j = 0, len = rest.length; j < len; j++) {
        a = rest[j];
        previous = this.derived_blade(previous, a[0], a[1], a[2]);
        if (!previous) {
          break;
        }
        results.push(previous);
      }
      return results;
    }).call(this);
    if (rest.length !== (config.length - 1)) {
      return null;
    }
    return [first].concat(rest);
  }

  combination_tests() {
    var a, combinations, config1, config2, g1, g2, grades1, grades2, j, json, l, len, len1, len2, len3, len4, len5, len6, len7, m, mv, n, o, p, q, r, ref, ref1, ref2, results, s, t1, t2;
    combinations = [];
    ref = [[0, "s"], [1, "e"]];
    for (j = 0, len = ref.length; j < len; j++) {
      [g1, t1] = ref[j];
      ref1 = [[0, "s"], [1, "e"]];
      for (l = 0, len1 = ref1.length; l < len1; l++) {
        [g2, t2] = ref1[l];
        if ("e" === t1 && t2 === t1) {
          continue;
        }
        combinations.push(this.mv([[g1, t1, 0], [g2, t2, 0]]));
      }
    }
    config1 = [[[1], "e"], [[1, 2], "n"], [[2], "en"]];
    config2 = [[[1, 2], "e"], [[1, 2], "n"], [[2, 3], "en"]];
    for (m = 0, len2 = config1.length; m < len2; m++) {
      [grades1, t1] = config1[m];
      for (n = 0, len3 = grades1.length; n < len3; n++) {
        g1 = grades1[n];
        for (p = 0, len4 = config2.length; p < len4; p++) {
          [grades2, t2] = config2[p];
          for (q = 0, len5 = grades2.length; q < len5; q++) {
            g2 = grades2[q];
            ref2 = [0, 2];
            for (r = 0, len6 = ref2.length; r < len6; r++) {
              o = ref2[r];
              mv = this.mv([[g1, t1, 0], [g2, t2, o]]);
              if (!mv) {
                continue;
              }
              combinations.push(mv);
            }
          }
        }
      }
    }
    results = [];
    for (s = 0, len7 = combinations.length; s < len7; s++) {
      a = combinations[s];
      json = JSON.stringify([a, this.space.mv_string(this.space.mv(a))]);
      results.push(console.log(json.replace(/,/g, ", ")));
    }
    return results;
  }

};

c3 = new sph_ga([1, 1, 1], {
  conformal: true
});

console.log(c3);

runner = new test_runner_class(c3);

gen = new test_data_generator_class(c3);

gen.combination_tests();

process.exit(0);

basic_tests = [
  {
    title: "id_bit_indices",
    actual: function() {
      var i,
  j,
  ref,
  results;
      results = [];
      for (i = j = 0, ref = 2 ** c3.n; (0 <= ref ? j < ref : j > ref); i = 0 <= ref ? ++j : --j) {
        results.push(c3.id_bit_indices(i));
      }
      return results;
    },
    expected: [[],
  [0],
  [1],
  [0,
  1],
  [2],
  [0,
  2],
  [1,
  2],
  [0,
  1,
  2]]
  }
];

map_grade_factor_tests = [
  {
    title: "Reverse: reverse(1) == 1",
    actual: function() {
      return c3.reverse(c3.vector([1]));
    },
    expected: c3.vector([1])
  },
  {
    title: "Reverse: reverse(e1) == e1",
    actual: function() {
      return c3.reverse(e1);
    },
    expected: e1
  },
  {
    title: "Reverse: reverse(e1e2) == -e1e2",
    actual: function() {
      return c3.reverse(c3.mv([[[1,
  2],
  1]]));
    },
    expected: c3.mv([[[1,
  2],
  -1]])
  },
  {
    title: "Reverse: reverse(e1e2e3) == -e1e2e3",
    actual: function() {
      return c3.reverse(c3.mv([[[1,
  2,
  3],
  1]]));
    },
    expected: c3.mv([[[1,
  2,
  3],
  -1]])
  },
  {
    title: "Reverse: reverse(1 + e1) == 1 + e1",
    actual: function() {
      return c3.reverse(c3.vector([1,
  1]));
    },
    expected: c3.vector([1,
  1])
  },
  {
    title: "Reverse: reverse(1 + e1 + e1e2) == 1 + e1 - e1e2",
    actual: function() {
      return c3.reverse(c3.mv([[[0],
  1],
  [[1],
  1],
  [[1,
  2],
  1]]));
    },
    expected: c3.mv([[[0],
  1],
  [[1],
  1],
  [[1,
  2],
  -1]])
  },
  {
    title: "Reverse: reverse(e1 + e2 + e1e2) == e1 + e2 - e1e2",
    actual: function() {
      return c3.reverse(c3.mv([[[1],
  1],
  [[2],
  1],
  [[1,
  2],
  1]]));
    },
    expected: c3.mv([[[1],
  1],
  [[2],
  1],
  [[1,
  2],
  -1]])
  },
  {
    title: "Reverse: reverse(1 + e1 + e2 + e1e2 + e1e3 + e1e2e3) == 1 + e1 + e2 - e1e2 - e1e3 - e1e2e3",
    actual: function() {
      return c3.reverse(c3.mv([
        [
          [0],
          1 // Scalar: 1
        ],
        [
          [1],
          1 // e1
        ],
        [
          [2],
          1 // e2
        ],
        [
          [1,
          2],
          1 // e1e2
        ],
        [
          [1,
          3],
          1 // e1e3
        ],
        [
          [1,
          2,
          3],
          1 // e1e2e3
        ]
      ]));
    },
    expected: c3.mv([
      [
        [0],
        1 // Scalar: 1
      ],
      [
        [1],
        1 // e1
      ],
      [
        [2],
        1 // e2
      ],
      [
        [1,
        2],
        -1 // -e1e2
      ],
      [
        [1,
        3],
        -1 // -e1e3
      ],
      [
        [1,
        2,
        3],
        -1 // -e1e2e3
      ]
    ])
  },
  {
    title: "Reverse Twice: reverse(reverse(e1 + e1e2)) == e1 + e1e2",
    actual: function() {
      return c3.reverse(c3.reverse(c3.mv([[[1],
  1],
  [[1,
  2],
  1]])));
    },
    expected: c3.mv([[[1],
  1],
  [[1,
  2],
  1]])
  },
  {
    title: "Involute: involute(e1e2e3) == -e1e2e3",
    actual: function() {
      return c3.involute(c3.mv([[[1,
  2,
  3],
  1]]));
    },
    expected: c3.mv([[[1,
  2,
  3],
  -1]])
  },
  {
    title: "Conjugate: conjugate(e1e2) == -e1e2",
    actual: function() {
      return c3.conjugate(c3.mv([[[1,
  2],
  1]]));
    },
    expected: c3.mv([[[1,
  2],
  -1]])
  }
];

ep_tests = [
  {
    title: "ep: Anticommutativity: a ∧ b == -b ∧ a",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          0,
          0 // a = e1
        ]),
        b: c3.vector([
          0,
          0,
          1,
          0 // b = e2
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.negate(c3.ep(d.b,
  d.a));
    }
  },
  {
    title: "ep: Exterior product of vector with itself is zero: a ∧ a == 0",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          2,
          3 // a = e1 + 2e2 + 3e3
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.a);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "ep: Associativity: (a ∧ b) ∧ c == a ∧ (b ∧ c)",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          0,
          0 // a = e1
        ]),
        b: c3.vector([
          0,
          0,
          1,
          0 // b = e2
        ]),
        c: c3.vector([
          0,
          0,
          0,
          1 // c = e3
        ])
      };
    },
    actual: function(d) {
      return c3.ep(c3.ep(d.a,
  d.b),
  d.c);
    },
    expected: function(d) {
      return c3.ep(d.a,
  c3.ep(d.b,
  d.c));
    }
  },
  {
    title: "ep: Distributivity over addition: a ∧ (b + c) == a ∧ b + a ∧ c",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          0,
          0 // a = e1
        ]),
        b: c3.vector([
          0,
          0,
          1,
          0 // b = e2
        ]),
        c: c3.vector([
          0,
          0,
          0,
          1 // c = e3
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  c3.add(d.b,
  d.c));
    },
    expected: function(d) {
      return c3.add(c3.ep(d.a,
  d.b),
  c3.ep(d.a,
  d.c));
    }
  },
  {
    title: "ep: Scalar multiplication: s ∧ a == s * a",
    data: function() {
      return {
        s: c3.s(2),
        a: c3.vector([
          0,
          1,
          0,
          0 // a = e1
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.s,
  d.a);
    },
    expected: function(d) {
      return c3.scale(d.a,
  2);
    }
  },
  {
    title: "ep: Exterior product of scalars: s ∧ t == 0",
    data: function() {
      return {
        s: c3.s(2),
        t: c3.s(3)
      };
    },
    actual: function(d) {
      return c3.ep(d.s,
  d.t);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "ep: Exterior product of bivectors: (e1 ∧ e2) ∧ (e2 ∧ e3) == 0",
    data: function() {
      return {
        b1: c3.mv([
          [
            [1,
            2],
            1 // b1 = e1 ∧ e2
          ]
        ]),
        b2: c3.mv([
          [
            [2,
            3],
            1 // b2 = e2 ∧ e3
          ]
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.b1,
  d.b2);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "ep: Exterior product of vector and bivector: a ∧ (e2 ∧ e3)",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          0,
          0 // a = e1
        ]),
        b: c3.mv([
          [
            [2,
            3],
            1 // b = e2 ∧ e3
          ]
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.mv([
        [
          [1,
          2,
          3],
          1 // Expected: e1 ∧ e2 ∧ e3
        ]
      ]);
    }
  },
  {
    title: "ep: Scalar times bivector: s ∧ (e1 ∧ e2) == s * (e1 ∧ e2)",
    data: function() {
      return {
        s: c3.s(2),
        b: c3.mv([
          [
            [1,
            2],
            1 // b = e1 ∧ e2
          ]
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.s,
  d.b);
    },
    expected: function(d) {
      return c3.scale(d.b,
  2);
    }
  },
  {
    title: "ep: Zero result when common basis elements: e1 ∧ e1 == 0",
    data: function() {
      return {
        e1: c3.basis(1)
      };
    },
    actual: function(d) {
      return c3.ep(d.e1,
  d.e1);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "cga: ep: Exterior product of two scalars is zero",
    data: function() {
      return {
        s1: c3.s(2),
        s2: c3.s(3)
      };
    },
    actual: function(d) {
      return c3.ep(d.s1,
  d.s2);
    },
    expected: function(d) {
      return c3.mv([]);
    }
  },
  {
    title: "cga: ep: Exterior product of a vector with itself is zero",
    data: function() {
      return {
        a: c3.vector([
          0,
          1,
          2,
          3,
          0,
          0 // a = 1e1 + 2e2 + 3e3
        ])
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.a);
    },
    expected: function(d) {
      return c3.mv([]);
    }
  },
  {
    title: "cga: ep: Exterior product of e1 and e2 equals e12",
    data: function() {
      return {
        e1: c3.basis(1),
        e2: c3.basis(2)
      };
    },
    actual: function(d) {
      return c3.ep(d.e1,
  d.e2);
    },
    expected: function(d) {
      return c3.mv([[[1,
  2],
  1]]);
    }
  },
  {
    title: "cga: ep: Antisymmetry: a ∧ b == - (b ∧ a)",
    data: function() {
      return {
        a: c3.basis(1),
        b: c3.basis(2)
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.scale(c3.ep(d.b,
  d.a),
  -1);
    }
  },
  {
    title: "cga: ep: Associativity: a ∧ (b ∧ c) == (a ∧ b) ∧ c",
    data: function() {
      return {
        a: c3.basis(1),
        b: c3.basis(2),
        c: c3.basis(3)
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  c3.ep(d.b,
  d.c));
    },
    expected: function(d) {
      return c3.ep(c3.ep(d.a,
  d.b),
  d.c);
    }
  },
  {
    title: "cga: ep: Exterior product of eo and ei equals eo ∧ ei",
    data: function() {
      return {
        eo: c3.eo(1),
        ei: c3.ei(1)
      };
    },
    actual: function(d) {
      return c3.ep(d.eo,
  d.ei);
    },
    expected: function(d) {
      return c3.mv([[[4,
  5],
  1]]);
    }
  },
  {
    title: "cga: ep: Exterior product of eo with itself is zero",
    data: function() {
      return {
        eo: c3.basis(4)
      };
    },
    actual: function(d) {
      return c3.ep(d.eo,
  d.eo);
    },
    expected: function(d) {
      return c3.mv([]);
    }
  },
  {
    title: "cga: ep: Exterior product of ei with itself is zero",
    data: function() {
      return {
        ei: c3.basis(5)
      };
    },
    actual: function(d) {
      return c3.ep(d.ei,
  d.ei);
    },
    expected: function(d) {
      return c3.mv([]);
    }
  },
  {
    title: "cga: ep: Exterior product of eo and e1 equals eo ∧ e1",
    data: function() {
      return {
        eo: c3.basis(4),
        e1: c3.basis(1)
      };
    },
    actual: function(d) {
      return c3.ep(d.eo,
  d.e1);
    },
    expected: function(d) {
      return c3.mv([[[1,
  4],
  -1]]);
    }
  },
  {
    title: "cga: ep: Antisymmetry with null vector: e1 ∧ eo == - (eo ∧ e1)",
    data: function() {
      return {
        e1: c3.basis(1),
        eo: c3.eo(4)
      };
    },
    actual: function(d) {
      return c3.ep(d.e1,
  d.eo);
    },
    expected: function(d) {
      return c3.mv([[[1,
  c3.eo_index],
  4]]);
    }
  },
  {
    title: "cga: ep: Exterior product of a scalar and a vector equals scalar times vector",
    data: function() {
      return {
        s1: c3.s(2),
        v: c3.vector([0,
  1,
  2,
  3,
  0,
  0])
      };
    },
    actual: function(d) {
      return c3.ep(d.s1,
  d.v);
    },
    expected: function(d) {
      return c3.scale(d.v,
  2);
    }
  },
  {
    title: "cga: ep: Exterior product of a vector and a scalar equals scalar times vector",
    data: function() {
      return {
        s1: c3.s(2),
        v: c3.vector([0,
  1,
  2,
  3])
      };
    },
    actual: function(d) {
      return c3.ep(d.v,
  d.s1);
    },
    expected: function(d) {
      return c3.scale(d.v,
  2);
    }
  },
  {
    title: "cga: ep: Exterior product of a blade with itself is zero",
    data: function() {
      return {
        a: c3.ep(c3.basis(1),
  c3.basis(2))
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.a);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "cga: ep: Exterior product of linearly dependent vectors is zero",
    data: function() {
      return {
        a: c3.basis(1),
        b: c3.basis(1,
  2)
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.s(0);
    }
  },
  {
    title: "cga: ep: Exterior product of a vector and a bivector",
    data: function() {
      return {
        a: c3.basis(1),
        b: c3.ep(c3.basis(2),
  c3.basis(3))
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.mv([[[1,
  2,
  3],
  1]]);
    }
  },
  {
    title: "cga: ep: Antisymmetry: a ∧ b == - (b ∧ a) for vector and bivector",
    data: function() {
      return {
        a: c3.basis(1),
        b: c3.ep(c3.basis(3),
  c3.basis(2))
      };
    },
    actual: function(d) {
      return c3.ep(d.a,
  d.b);
    },
    expected: function(d) {
      return c3.mv([[[1,
  2,
  3],
  -1]]);
    }
  }
];

ip_tests = [
  {
    title: "ip: scalar ⋅ scalar",
    actual: function() {
      return c3.ip(s2,
  s3);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: scalar ⋅ mv_nn",
    actual: function() {
      return c3.ip(s2,
  e1);
    },
    expected: c3.basis(1,
  4)
  },
  {
    title: "ip: mv_nn ⋅ scalar",
    actual: function() {
      return c3.ip(e1,
  s3);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: scalar ⋅ mv_n",
    actual: function() {
      return c3.ip(s2,
  eo);
    },
    expected: c3.eo(10)
  },
  {
    title: "ip: mv_n ⋅ scalar",
    actual: function() {
      return c3.ip(eo,
  s3);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_n1 ⋅ mv_n2",
    actual: function() {
      return c3.ip(eo,
  ei);
    },
    expected: c3.s(-30)
  },
  {
    title: "ip: mv_n1 ⋅ mv_n1",
    actual: function() {
      return c3.ip(eo,
  eo);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_n ⋅ mv_nn",
    actual: function() {
      return c3.ip(eo,
  e2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_nn ⋅ mv_n",
    actual: function() {
      return c3.ip(e2,
  eo);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_nn ⋅ mv_nn",
    actual: function() {
      return c3.ip(e1,
  e2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_nn ⋅ mv_nnn",
    actual: function() {
      return c3.ip(e1,
  vector_nnn);
    },
    expected: c3.s(4)
  },
  {
    title: "ip: mv_nnn ⋅ mv_nn",
    actual: function() {
      return c3.ip(vector_nnn,
  e2);
    },
    expected: c3.s(9)
  },
  {
    title: "ip: mv_nnn ⋅ mv_nnn",
    actual: function() {
      return c3.ip(vector_nnn,
  vector_nnn_ei);
    },
    expected: c3.s(-1)
  },
  {
    title: "ip: mv_n ⋅ mv_nnn",
    actual: function() {
      return c3.ip(eo,
  vector_nn);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_nnn ⋅ mv_n",
    actual: function() {
      return c3.ip(vector_nn,
  ei);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: bivector ⋅ vector",
    actual: function() {
      return c3.ip(bivector,
  e3);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: vector ⋅ bivector",
    actual: function() {
      return c3.ip(e3,
  bivector);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: trivector ⋅ bivector",
    actual: function() {
      return c3.ip(trivector,
  bivector);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: mv_n ⋅ bivector",
    actual: function() {
      return c3.ip(eo,
  bivector);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: bivector ⋅ mv_n",
    actual: function() {
      return c3.ip(bivector,
  ei);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e1 ⋅ e1", // Full basis overlap
    actual: function() {
      return c3.ip(e1,
  e1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e1 ⋅ e2", // No basis overlap
    actual: function() {
      return c3.ip(e1,
  e2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e1 ⋅ vector_nn1", // Partial basis overlap
    actual: function() {
      return c3.ip(e1,
  vector_nn1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: vector_nn1 ⋅ vector_nn1", // Full basis overlap
    actual: function() {
      return c3.ip(vector_nn1,
  vector_nn1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: vector_nn1 ⋅ vector_nn2", // Partial basis overlap
    actual: function() {
      return c3.ip(vector_nn1,
  vector_nn2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: bivector1 ⋅ bivector1", // Full basis overlap
    actual: function() {
      return c3.ip(bivector1,
  bivector1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: bivector1 ⋅ bivector2", // Partial basis overlap
    actual: function() {
      return c3.ip(bivector1,
  bivector2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: bivector1 ⋅ bivector3", // No basis overlap
    actual: function() {
      return c3.ip(bivector1,
  bivector3);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e1 ⋅ bivector1", // Partial basis overlap
    actual: function() {
      return c3.ip(e1,
  bivector1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e1 ⋅ bivector2", // Partial basis overlap
    actual: function() {
      return c3.ip(e1,
  bivector2);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: e2 ⋅ bivector1", // Partial basis overlap
    actual: function() {
      return c3.ip(e2,
  bivector1);
    },
    expected: c3.s(0)
  },
  {
    title: "ip: trivector1 ⋅ trivector1", // Full basis overlap
    actual: function() {
      return c3.ip(trivector1,
  trivector1);
    },
    expected: c3.s(0)
  }
];

runner.run_tests(ip_tests);

//basic_tests
//map_grade_factor_tests
//ep_tests
