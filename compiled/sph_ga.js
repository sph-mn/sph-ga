// Generated by CoffeeScript 2.7.0
var sph_ga,
  indexOf = [].indexOf;

sph_ga = (function() {
  class sph_ga {
    constructor(metric, options = {}) {
      var ref, ref1, rotation_axis_combinations;
      this.n = metric.length;
      this.is_conformal = !!options.conformal;
      if (!Array.isArray(metric[0])) {
        if (this.is_conformal) {
          this.n += 2;
          metric.push(0, 0);
        }
        metric = this.flat_metric_to_full(metric, this.n);
        if (this.is_conformal) {
          metric[this.n - 1][this.n - 2] = -1;
          metric[this.n - 2][this.n - 1] = -1;
        }
      }
      [this.is_symmetric, this.is_diagonal, this.null_vector_count] = this.metric_properties(metric, this.n);
      if (!this.is_symmetric) {
        throw new Error("the metric must be symmetric in the non-null part");
      }
      if (this.is_conformal) {
        if (2 !== this.null_vector_count) {
          throw new Error("only two null vectors are allowed with \"conformal: true\". use a custom metric instead");
        }
        if (!this.is_diagonal) {
          throw new Error("only diagonal metrics are allowed with \"conformal: true\". use a custom metric instead");
        }
      }
      if (this.null_vector_count) {
        this.null_vector_start = this.n - this.null_vector_count;
        this.id_null = this.id_from_indices((function() {
          var results = [];
          for (var l = ref = this.null_vector_start + 1, ref1 = this.n; ref <= ref1 ? l <= ref1 : l >= ref1; ref <= ref1 ? l++ : l--){ results.push(l); }
          return results;
        }).apply(this));
      }
      this.pseudoscalar_id = (1 << this.n) - 1;
      this.metric = metric;
      if (this.is_diagonal) {
        this.ip_metric = function(indices) {
          var i;
          return this.array_product((function() {
            var l, len, results;
            results = [];
            for (l = 0, len = indices.length; l < len; l++) {
              i = indices[l];
              results.push(this.metric[i][i]);
            }
            return results;
          }).call(this));
        };
      } else {
        this.ip_metric = function(indices) {
          var i, j;
          if (!indices.length) {
            return 1;
          }
          return this.determinant((function() {
            var l, len, results;
            results = [];
            for (l = 0, len = indices.length; l < len; l++) {
              i = indices[l];
              results.push((function() {
                var len1, o, results1;
                results1 = [];
                for (o = 0, len1 = indices.length; o < len1; o++) {
                  j = indices[o];
                  results1.push(this.metric[i][j]);
                }
                return results1;
              }).call(this));
            }
            return results;
          }).call(this));
        };
      }
      if (this.is_conformal) {
        this.no_bit_index = this.n - 2;
        this.ni_bit_index = this.n - 1;
        this.no_id = 1 << this.no_bit_index;
        this.ni_id = 1 << this.ni_bit_index;
        this.no_index = this.no_bit_index + 1;
        this.ni_index = this.ni_bit_index + 1;
        this.no = function(coeff) {
          return [[this.no_id, coeff, 1]];
        };
        this.ni = function(coeff) {
          return [[this.ni_id, coeff, 1]];
        };
        rotation_axis_combinations = function(n) {
          var combinations, i, j, l, o, ref2, ref3, ref4;
          combinations = [];
          for (i = l = 0, ref2 = n; (0 <= ref2 ? l < ref2 : l > ref2); i = 0 <= ref2 ? ++l : --l) {
            for (j = o = ref3 = i + 1, ref4 = n; (ref3 <= ref4 ? o < ref4 : o > ref4); j = ref3 <= ref4 ? ++o : --o) {
              combinations.push([i + 1, j + 1]);
            }
          }
          return combinations;
        };
        this.rotation_axes = rotation_axis_combinations(this.n - 2);
        this.rotor = function(coeffs) {
          var a, first, i, rest;
          [first, ...rest] = coeffs;
          return this.mv([[[0], first]].concat((function() {
            var l, len, results;
            results = [];
            for (i = l = 0, len = rest.length; l < len; i = ++l) {
              a = rest[i];
              results.push([this.rotation_axes[i], a]);
            }
            return results;
          }).call(this)));
        };
        this.point = function(euclidean_coeffs) {
          var a, ni_coeff;
          ni_coeff = this.array_sum((function() {
            var l, len, results;
            results = [];
            for (l = 0, len = euclidean_coeffs.length; l < len; l++) {
              a = euclidean_coeffs[l];
              results.push(a * a);
            }
            return results;
          })()) / 2;
          return this.vector([0].concat(euclidean_coeffs).concat([1, ni_coeff]));
        };
        this.point_euclidean = function(a) {
          var b, c, d, l, len, n0, results;
          n0 = this.blade_coeff(this.get(a, this.no_id));
          c = (function() {
            var l, len, results;
            results = [];
            for (l = 0, len = a.length; l < len; l++) {
              b = a[l];
              if (1 === b[2] && !(b[0] === this.no_id || b[0] === this.ni_id)) {
                results.push(b[1]);
              }
            }
            return results;
          }).call(this);
          results = [];
          for (l = 0, len = c.length; l < len; l++) {
            d = c[l];
            results.push(d / n0);
          }
          return results;
        };
        this.normal = this.vector([0].concat(Array(this.n).fill(1 / Math.sqrt(this.n))));
      }
    }

    add_one(a, b) {
      return this.combine(a, b, 1);
    }

    array_diff(a, b) {
      return a.filter(function(c) {
        return !(indexOf.call(b, c) >= 0);
      });
    }

    array_product(a) {
      return a.reduce((function(b, a) {
        return a * b;
      }), 1);
    }

    array_sum(a) {
      return a.reduce((function(b, a) {
        return a + b;
      }), 0);
    }

    basis_blade(i, coeff) {
      if (i) {
        return [1 << (i - 1), coeff, 1];
      } else {
        return [0, coeff, 0];
      }
    }

    basis(i, coeff = 1) {
      return [this.basis_blade(i, coeff)];
    }

    blade_coeff(a) {
      return a[1];
    }

    blade_grade(a) {
      return a[2];
    }

    blade_id(a) {
      return a[0];
    }

    coeffs_add(coeffs, id, coeff, grade) {
      if (coeffs[id] != null) {
        return coeffs[id][0] += coeff;
      } else {
        return coeffs[id] = [coeff, grade];
      }
    }

    conjugate(a) {
      return this.map_grade_factor(a, function(grade) {
        return (-1) ** ((grade * (grade + 1)) >> 1);
      });
    }

    get(a, id) {
      var b, l, len;
      for (l = 0, len = a.length; l < len; l++) {
        b = a[l];
        if (id === b[0]) {
          return b;
        }
      }
    }

    grade(a) {
      return a[a.length - 1][2];
    }

    id_from_bit_indices(indices) {
      return indices.reduce((function(id, i) {
        return id |= 1 << i;
      }), 0);
    }

    id_from_indices(indices) {
      return indices.reduce((function(id, i) {
        return id |= 1 << (i - 1);
      }), 0);
    }

    involute(a) {
      return this.map_grade_factor(a, function(grade) {
        return (-1) ** grade;
      });
    }

    map_grade_factor(a, f) {
      var coeff, grade, id, l, len, results;
      results = [];
      for (l = 0, len = a.length; l < len; l++) {
        [id, coeff, grade] = a[l];
        results.push([id, coeff * f(grade), grade]);
      }
      return results;
    }

    mv(terms) {
      var coeff, indices, l, len, results;
      results = [];
      for (l = 0, len = terms.length; l < len; l++) {
        [indices, coeff] = terms[l];
        results.push(this.blade(indices, coeff));
      }
      return results;
    }

    negate(a) {
      return this.scale(a, -1);
    }

    pseudoscalar() {
      var ref;
      return [
        this.blade((function() {
          var results = [];
          for (var l = 1, ref = this.n; 1 <= ref ? l <= ref : l >= ref; 1 <= ref ? l++ : l--){ results.push(l); }
          return results;
        }).apply(this),
        1)
      ];
    }

    reverse(a) {
      return this.map_grade_factor(a, function(grade) {
        return (-1) ** ((grade * (grade - 1)) >> 1);
      });
    }

    scale(mv, a) {
      var coeff, grade, id, l, len, results;
      results = [];
      for (l = 0, len = mv.length; l < len; l++) {
        [id, coeff, grade] = mv[l];
        results.push([id, coeff * a, grade]);
      }
      return results;
    }

    s(coeff) {
      return [[0, coeff, 0]];
    }

    sp(a, b) {
      return this.gp(this.gp(a, b), this.inverse(a));
    }

    subtract_one(a, b) {
      return this.combine(a, b, -1);
    }

    vector(coeffs) {
      var a, i, l, len, results;
      results = [];
      for (i = l = 0, len = coeffs.length; l < len; i = ++l) {
        a = coeffs[i];
        if (a) {
          results.push(this.basis_blade(i, a));
        }
      }
      return results;
    }

    ep(...a) {
      return a.reduce((c, b) => {
        return this.ep_one(c, b);
      });
    }

    ip(...a) {
      return a.reduce((c, b) => {
        return this.ip_one(c, b);
      });
    }

    gp(...a) {
      return a.reduce((c, b) => {
        return this.gp_one(c, b);
      });
    }

    add(...a) {
      return a.reduce((c, b) => {
        return this.add_one(c, b);
      });
    }

    subtract(...a) {
      return a.reduce((c, b) => {
        return this.subtract_one(c, b);
      });
    }

    mv_to_string(a) {
      var b;
      return ((function() {
        var l, len, results;
        results = [];
        for (l = 0, len = a.length; l < len; l++) {
          b = a[l];
          results.push(this.blade_to_string(b));
        }
        return results;
      }).call(this)).join(" + ");
    }

    blade(indices, coeff) {
      if (indices[0]) {
        return [this.id_from_indices(indices), coeff, indices.length];
      } else {
        return [this.id_from_indices(indices.slice(1)), coeff, indices.length - 1];
      }
    }

    coeffs_to_mv(coeffs) {
      var a, coeff, grade, id;
      a = (function() {
        var results;
        results = [];
        for (id in coeffs) {
          [coeff, grade] = coeffs[id];
          if (coeff !== 0) {
            results.push([parseInt(id), coeff, grade]);
          }
        }
        return results;
      })();
      if (a.length) {
        return a;
      } else {
        return [[0, 0, 0]];
      }
    }

    id_indices(id) {
      var a, l, len, ref, results;
      if (id) {
        ref = this.id_bit_indices(id);
        results = [];
        for (l = 0, len = ref.length; l < len; l++) {
          a = ref[l];
          results.push(1 + a);
        }
        return results;
      } else {
        return [0];
      }
    }

    id_bit_indices(id) {
      var a, i;
      if (id in this.id_bit_indices_cache) {
        return this.id_bit_indices_cache[id];
      }
      a = (function() {
        var l, ref, results;
        results = [];
        for (i = l = 0, ref = this.n; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
          if (id & (1 << i)) {
            results.push(i);
          }
        }
        return results;
      }).call(this);
      this.id_bit_indices_cache[id] = a;
      return a;
    }

    flat_metric_to_full(metric, n) {
      var a, b, i, l, ref;
      a = Array(n);
      for (i = l = 0, ref = n; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
        b = Array(n).fill(0);
        b[i] = metric[i];
        a[i] = b;
      }
      return a;
    }

    metric_properties(metric, n) {
      var di, first_zero, i, is_diagonal, j, k, l, null_vectors, o, p, ref, ref1, ref2, ref3, ref4, ref5;
      null_vectors = 0;
      first_zero = n;
      is_diagonal = true;
      for (i = l = 0, ref = n; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
        di = metric[i][i];
        if (di === 0) {
          null_vectors += 1;
          if (first_zero === n) {
            first_zero = i;
          }
        } else {
          if (first_zero < n || ((ref1 = !di) === 1 || ref1 === (-1))) {
            is_diagonal = false;
            break;
          }
        }
        for (j = o = ref2 = i + 1, ref3 = n; (ref2 <= ref3 ? o < ref3 : o > ref3); j = ref2 <= ref3 ? ++o : --o) {
          if (i < first_zero) {
            if (metric[i][j] !== metric[j][i]) {
              return [false, false, null_vectors];
            }
            if (metric[i][j] !== 0) {
              is_diagonal = false;
              break;
            }
          }
        }
        if (!is_diagonal) {
          break;
        }
      }
      for (k = p = ref4 = first_zero, ref5 = n; (ref4 <= ref5 ? p < ref5 : p > ref5); k = ref4 <= ref5 ? ++p : --p) {
        if (metric[k][k] !== 0) {
          return [false, false, null_vectors];
        }
      }
      return [true, is_diagonal, null_vectors];
    }

    id_grade(a) {
      var b, n;
      if (a in this.id_grade_cache) {
        return this.id_grade_cache[a];
      }
      n = 0;
      b = a;
      while (b !== 0) {
        b &= b - 1;
        n += 1;
      }
      this.id_grade_cache[a] = n;
      return n;
    }

    determinant_generic(a, n) {
      var b, c, i, j, l, o, ref, ref1, sign;
      if (n === 1) {
        return a[0][0];
      }
      b = 0;
      for (j = l = 0, ref = n; (0 <= ref ? l < ref : l > ref); j = 0 <= ref ? ++l : --l) {
        for (i = o = 1, ref1 = n; (1 <= ref1 ? o < ref1 : o > ref1); i = 1 <= ref1 ? ++o : --o) {
          c = a[i].slice(0, j).concat(a[i].slice(j + 1));
        }
        sign = 0 === j % 2 ? 1 : -1;
        b += sign * a[0][j] * this.determinant_generic(c);
      }
      return b;
    }

    determinant(matrix) {
      var n;
      n = matrix.length;
      if (1 === n) {
        return matrix[0][0];
      } else if (2 === n) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
      } else if (3 === n) {
        return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
      } else if (4 === n) {
        return matrix[0][0] * (matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) - matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) + matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])) - matrix[0][1] * (matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) - matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) + matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])) + matrix[0][2] * (matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) - matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) + matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])) - matrix[0][3] * (matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) - matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) + matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]));
      } else {
        return this.determinant_generic(matrix, n);
      }
    }

    array_remove_pair(a, i, j) {
      if (i < j) {
        a.splice(j, 1);
        return a.splice(i, 1);
      } else {
        a.splice(i, 1);
        return a.splice(j, 1);
      }
    }

    sign(indices) {
      var c, i, j, l, o, ref, ref1, ref2;
      c = 0;
      for (i = l = 0, ref = indices.length; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
        for (j = o = ref1 = i + 1, ref2 = indices.length; (ref1 <= ref2 ? o < ref2 : o > ref2); j = ref1 <= ref2 ? ++o : --o) {
          if (indices[i] > indices[j]) {
            c += 1;
          }
        }
      }
      return (-1) ** c;
    }

    sign_sorted(a, b) {
      var c, i, j;
      c = 0;
      i = 0;
      j = 0;
      while (i < a.length && j < b.length) {
        if (a[i] <= b[j]) {
          i += 1;
        } else {
          c += 1;
          j += 1;
        }
      }
      return (-1) ** c;
    }

    ip_one(a, b) {
      var c, choose, coeff_a, coeff_b, coeff_c, coeffs, grade_a, grade_b, grade_c, id_a, id_b, id_c, indices_a, indices_b, inner_product_terms, ip_equal_grade, l, len, len1, len2, o, p, sign_of_choice;
      choose = function(a, k) { // return all combinations of "k" elements from "a".
        var head, tail, with_head, without_head;
        if (k === 0) {
          return [[]];
        }
        if (a.length < k) {
          return [];
        }
        head = a[0];
        tail = a.slice(1);
        with_head = choose(tail, k - 1);
        with_head = with_head.map(function(combo) {
          return [head].concat(combo);
        });
        without_head = choose(tail, k);
        return with_head.concat(without_head);
      };
      sign_of_choice = function(combo) {
        var i, l, ref, s;
        s = 0;
        for (i = l = 0, ref = combo.length; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
          s += combo[i] - i;
        }
        return (-1) ** s;
      };
      inner_product_terms = (a, b) => {
        var c, chosen, coeff, combo, i, id, indices, l, len, n, o, ref, ref1, ref2, results;
        n = a.length;
        indices = (function() {
          var results = [];
          for (var l = 0, ref = b.length; 0 <= ref ? l < ref : l > ref; 0 <= ref ? l++ : l--){ results.push(l); }
          return results;
        }).apply(this);
        ref1 = choose(indices, n);
        results = [];
        for (l = 0, len = ref1.length; l < len; l++) {
          combo = ref1[l];
          chosen = (function() {
            var len1, o, results1;
            results1 = [];
            for (o = 0, len1 = combo.length; o < len1; o++) {
              i = combo[o];
              results1.push(b[i]);
            }
            return results1;
          })();
          c = indices.filter(function(i) {
            return !(indexOf.call(combo, i) >= 0);
          });
          coeff = sign_of_choice(combo);
          for (i = o = 0, ref2 = a.length; (0 <= ref2 ? o < ref2 : o > ref2); i = 0 <= ref2 ? ++o : --o) {
            coeff *= this.metric[a[i]][chosen[i]];
          }
          id = this.id_from_bit_indices((function() {
            var len1, p, results1;
            results1 = [];
            for (p = 0, len1 = c.length; p < len1; p++) {
              i = c[p];
              results1.push(b[i]);
            }
            return results1;
          })());
          results.push([id, coeff, c.length]);
        }
        return results;
      };
      ip_equal_grade = (a, b) => {
        var i, j, m, s;
        // for 1-blades use the metric directly.
        if (a.length === 1) {
          return this.s(this.metric[a[0]][b[0]]);
        } else {
          m = (function() {
            var l, ref, results;
// -> gram_matrix
            results = [];
            for (i = l = 0, ref = a.length; (0 <= ref ? l < ref : l > ref); i = 0 <= ref ? ++l : --l) {
              results.push((function() {
                var o, ref1, results1;
                results1 = [];
                for (j = o = 0, ref1 = b.length; (0 <= ref1 ? o < ref1 : o > ref1); j = 0 <= ref1 ? ++o : --o) {
                  results1.push(this.metric[a[i]][b[j]]);
                }
                return results1;
              }).call(this));
            }
            return results;
          }).call(this);
          s = (-1) ** (a.length * (a.length - 1) / 2);
          return this.s(-s * this.determinant(m));
        }
      };
      coeffs = {};
      for (l = 0, len = a.length; l < len; l++) {
        [id_a, coeff_a, grade_a] = a[l];
        indices_a = this.id_bit_indices(id_a);
        for (o = 0, len1 = b.length; o < len1; o++) {
          [id_b, coeff_b, grade_b] = b[o];
          if (grade_b < grade_a) {
            continue;
          }
          indices_b = this.id_bit_indices(id_b);
          if (grade_a === grade_b) {
            c = ip_equal_grade(indices_a, indices_b);
          } else {
            c = inner_product_terms(indices_a, indices_b);
          }
          for (p = 0, len2 = c.length; p < len2; p++) {
            [id_c, coeff_c, grade_c] = c[p];
            this.coeffs_add(coeffs, id_c, coeff_a * coeff_b * coeff_c, grade_c);
          }
        }
      }
      return this.coeffs_to_mv(coeffs);
    }

    gp_one(a, b) {
      var coeff_a, coeff_b, coeffs, combined, final_coeff, grade_a, grade_b, id_a, id_b, indices_a, indices_b, l, len, len1, m_factor, merge_indices, merged, o, sign;
      merge_indices = (indices) => {
        var changed, factor, i, j, m;
        // merge the basis vector indices from two blades while incorporating metric factors.
        factor = 1;
        changed = true;
        // first pass: process pairs of distinct indices with non-zero off-diagonal metric terms.
        while (changed) {
          changed = false;
          i = 0;
          while (i < indices.length) {
            j = i + 1;
            while (j < indices.length) {
              // check if indices differ; if so, see if their off-diagonal metric element contributes.
              if (indices[i] !== indices[j]) {
                m = this.metric[indices[i]][indices[j]];
                if (m !== 0) {
                  factor *= m;
                  // remove the pair because their metric factor has been factored out.
                  this.array_remove_pair(indices, i, j);
                  changed = true;
                  break;
                } else {
                  j += 1;
                }
              } else {
                j += 1;
              }
            }
            i += 1;
          }
        }
        changed = true;
        // second pass: process duplicate indices using the corresponding diagonal metric factors.
        while (changed) {
          changed = false;
          i = 0;
          while (i < indices.length) {
            j = i + 1;
            while (j < indices.length) {
              // when two identical indices are found, multiply by the square from the metric.
              if (indices[i] === indices[j]) {
                m = this.metric[indices[i]][indices[i]];
                factor *= m;
                // remove the duplicate pair to simplify the blade.
                this.array_remove_pair(indices, i, j);
                changed = true;
                break;
              } else {
                j += 1;
              }
            }
            i += 1;
          }
        }
        return [indices, factor];
      };
      coeffs = {};
      for (l = 0, len = a.length; l < len; l++) {
        [id_a, coeff_a, grade_a] = a[l];
        indices_a = this.id_bit_indices(id_a);
        for (o = 0, len1 = b.length; o < len1; o++) {
          [id_b, coeff_b, grade_b] = b[o];
          indices_b = this.id_bit_indices(id_b);
          combined = indices_a.concat(indices_b);
          sign = this.sign(combined);
          [merged, m_factor] = merge_indices(combined.slice());
          final_coeff = coeff_a * coeff_b * sign * m_factor;
          if (!final_coeff) {
            continue;
          }
          this.coeffs_add(coeffs, this.id_from_bit_indices(merged), final_coeff, merged.length);
        }
      }
      return this.coeffs_to_mv(coeffs);
    }

    ep_one(a, b) {
      var coeff_a, coeff_b, coeffs, grade_a, grade_b, id, id_a, id_b, indices_a, l, len, len1, o, sign;
      coeffs = {};
      for (l = 0, len = a.length; l < len; l++) {
        [id_a, coeff_a, grade_a] = a[l];
        indices_a = this.id_bit_indices(id_a);
        for (o = 0, len1 = b.length; o < len1; o++) {
          [id_b, coeff_b, grade_b] = b[o];
          if (!(!(id_a & id_b))) {
            continue;
          }
          id = id_a | id_b;
          if (!id) {
            continue;
          }
          sign = id_b ? this.sign_sorted(indices_a, this.id_bit_indices(id_b)) : 1;
          this.coeffs_add(coeffs, id, sign * coeff_a * coeff_b, grade_a + grade_b);
        }
      }
      return this.coeffs_to_mv(coeffs);
    }

    combine(a, b, scalar = 1) {
      var c, coeff, coeffs, grade, id, l, len, len1, o;
      coeffs = {};
      for (l = 0, len = a.length; l < len; l++) {
        [id, coeff, grade] = a[l];
        coeffs[id] = [coeff, grade];
      }
      for (o = 0, len1 = b.length; o < len1; o++) {
        [id, coeff, grade] = b[o];
        if (coeffs[id] != null) {
          coeffs[id][0] += coeff * scalar;
        } else {
          coeffs[id] = [coeff * scalar, grade];
        }
      }
      c = (function() {
        var results;
        results = [];
        for (id in coeffs) {
          [coeff, grade] = coeffs[id];
          if (coeff !== 0) {
            results.push([parseInt(id), coeff, grade]);
          }
        }
        return results;
      })();
      if (c.length) {
        return c;
      } else {
        return this.null_scalar;
      }
    }

    inverse(a) {
      var a_reverse, coeff, denom, denom_mv, grade, id, l, len, len1, o, results;
      a_reverse = this.reverse(a);
      denom_mv = this.gp(a, a_reverse);
      denom = 0;
      for (l = 0, len = denom_mv.length; l < len; l++) {
        [id, coeff, grade] = denom_mv[l];
        if (id === 0) {
          denom = coeff;
          break;
        }
      }
      if (denom === 0) {
        throw new Error("multivector is not invertible (denominator is zero).");
      }
      results = [];
      for (o = 0, len1 = a_reverse.length; o < len1; o++) {
        [id, coeff, grade] = a_reverse[o];
        results.push([parseInt(id), coeff / denom, grade]);
      }
      return results;
    }

    blade_to_string(a) {
      var base, coeff, grade, id;
      [id, coeff, grade] = a;
      if (id) {
        base = "e" + this.id_bit_indices(id).map(function(b) {
          return b + 1;
        }).join("_");
        if (1 === coeff) {
          return base;
        } else if (-1 === coeff) {
          return `-${base}`;
        } else {
          return coeff + base;
        }
      } else {
        return coeff;
      }
    }

    blade_from_string(a) {
      var coeff, id, indices, left_number, letters, match, n, result, right_numbers;
      match = a.match(/^(?:(\d+(?:\.\d+)?))?([a-z]+)?(?:([\d_]+))?$/);
      if (!match) {
        return null;
      }
      left_number = match[1] != null ? parseFloat(match[1]) : null;
      letters = match[2] != null ? match[2] : null;
      right_numbers = match[3] != null ? (function() {
        var l, len, ref, results;
        ref = match[3].split("_");
        results = [];
        for (l = 0, len = ref.length; l < len; l++) {
          n = ref[l];
          results.push(parseInt(n));
        }
        return results;
      })() : null;
      result = [];
      coeff = left_number != null ? left_number : 1;
      indices = right_numbers != null ? right_numbers : [];
      if (letters != null) {
        switch (letters) {
          case "no":
            id = c3.no_id;
            break;
          case "ni":
            id = c3.ni_id;
            break;
          default:
            id = this.id_from_indices(right_numbers);
        }
      } else {
        id = 0;
      }
      return [id, coeff, this.id_grade(id)];
    }

    mv_from_string(a) {
      var b, l, len, ref, results;
      ref = a.split(" + ");
      results = [];
      for (l = 0, len = ref.length; l < len; l++) {
        b = ref[l];
        results.push(this.blade_from_string(b));
      }
      return results;
    }

  };

  sph_ga.prototype.id_bit_indices_cache = {};

  sph_ga.prototype.id_grade_cache = {};

  sph_ga.prototype.null_scalar = [[0, 0, 0]];

  return sph_ga;

}).call(this);

if (typeof module !== "undefined" && module.exports) {
  module.exports = sph_ga;
} else {
  window.sph_ga = sph_ga;
}
