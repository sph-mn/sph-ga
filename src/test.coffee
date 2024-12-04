sph_ga = require "./ga.coffee"

class test_runner_class
  # this class executes automated tests based on an array format.
  is_string: (a) -> typeof a is "string"
  to_json: (a) -> JSON.stringify(a).replace /,(?=\S)/g, ", "
  is_plain_object: (a) -> a? and typeof a is "object" and a.constructor is Object
  any_to_array: (a) -> if Array.isArray(a) then a else [a]
  object_merge: (a, b) ->
    for k, v of b
      if @is_plain_object(v) and @is_plain_object(a[k])
        a[k] = @object_merge a[k], v
      else
        a[k] = v
    a
  report_compact_failure_strings: (inp, exp, out) -> [((@to_json a for a in inp).join ", "), @to_json(exp), @to_json(out)]
  report_compact: (results) ->
    for [name, test_results...] in results
      process.stdout.write name
      for [status, index, name, inp, exp, out] in test_results
        if status then process.stdout.write " #{index}"
        else
          [inp_string, exp_string, out_string] = @report_compact_failure_strings inp, exp, out
          process.stdout.write [
            "\n  failure #{name} #{index}"
            "inp #{inp_string}"
            "exp #{exp_string}"
            "out #{out_string}"
          ].join "\n  "
      console.log ""
  constructor: (options) ->
    default_options =
      reporter: @report_compact
    @options = @object_merge default_options, options
    @options.reporter = @options.reporter.bind @
  execute_tests: (tests) ->
    status = true
    for [f, rest...] in tests
      break unless status
      [name, context] =
        if Array.isArray f
          [name_or_context, f] = f
          if @is_string(name_or_context) then [name_or_context, null] else [f.name, name_or_context]
        else [f.name, null]
      results = [name]
      for i in [0...rest.length] by 2
        inp = @any_to_array rest[i]
        exp = rest[i + 1]
        out = f.apply context, inp
        out_string = @to_json out
        exp_string = @to_json exp
        status = out_string == exp_string
        results.push [status, i / 2, name, inp, exp, out]
        break unless status
      results
  execute: (tests) -> @options.reporter @execute_tests tests

class test_data_generator_class
  # this class systematically generates blade combinations and multivectors for testing.
  to_json: (a) -> JSON.stringify(a).replace /,(?=\S)/g, ", "
  random_integer: (min, max) -> Math.floor(Math.random() * (max - min + 1)) + min
  random_element: (a) -> a[Math.floor(Math.random() * a.length)]
  randomize: (a) ->
    for i in [(a.length - 1)..1]
      j = Math.floor(Math.random() * (i + 1))
      [a[i], a[j]] = [a[j], a[i]]
    a
  id_type: (id) ->
    if id & ~@space.id_null
      if id & @space.id_null then type = "en"
      else type = "e"
    else
      if id & @space.id_null then type = "n"
      else type = "s"
  constructor: (space) ->
    @cycled_coeff = @cycled_sequence [2..9]
    @space = space
    @n = space.n
    @blade_ids = ({s: [], e: [], n: [], en: []} for [0..@n])
    for id in [0...(2 ** @n)]
      g = space.id_grade id
      t = @id_type id
      @blade_ids[g][t].push id
    @blade_id_next = for a in @blade_ids
      b = {}
      for k, v of a
        b[k] = @cycled_sequence v
      b
  cycled_sequence: (a) ->
    i = 0
    ->
      b = a[i]
      i += 1
      i = if i < a.length then i else 0
      b
  blade: (grade, type, overlap_id) ->
    # type
    #   s: scalar
    #   e: euclidean
    #   n: null
    #   en: mixed
    #   any: any
    if "any" == type
      types = Object.keys(@blade_ids[grade]).filter (a) => @blade_ids[grade][a].length
      return null unless types.length
      type = @random_element types
    if overlap_id
      ids = @blade_ids[grade][type]
      ids_with_overlap = ids.filter (a) -> a & overlap_id
      id = @random_element ids_with_overlap
    else id = @blade_id_next[grade][type]()
    return null unless id?
    [@space.id_indices(id), @cycled_coeff()]
  derived_blade: (blade_a, grade_b, type_b, overlap_type) ->
    # overlap_type
    #   0: unspecified overlap
    #   1: partial overlap
    #   2: full overlap
    id_a = @space.id_from_indices blade_a[0]
    type_a = @id_type id_a
    grade_a = @space.id_grade id_a
    switch overlap_type
      when 0 then @blade grade_b, type_b
      when 1 then @blade grade_b, type_b, id_a
      when 2
        return null unless grade_a == grade_b && type_a == type_b
        blade_b = Array.from(blade_a)
        blade_b[1] = @cycled_coeff()
        blade_b
  mv: (config) ->
    # [[grade, type, overlap_type], ...]
    first = @blade.apply @, config[0][0..1]
    rest = config[1..]
    previous = first
    rest = for a in rest
      previous = @derived_blade previous, a[0], a[1], a[2]
      break unless previous
      previous
    return null if rest.length != (config.length - 1)
    [first].concat rest
  display_combination_tests: () ->
    combinations = []
    for [g1, t1] in [[0, "s"], [1, "e"]]
      for [g2, t2] in [[0, "s"], [1, "e"]]
        continue if "e" == t1 && t2 == t1
        mv = @mv [[g1, t1, 0], [g2, t2, 0]]
        combinations.push [[mv[0]], [mv[1]]]
    config1 = [[[1], "e"], [[1, 2], "n"], [[2], "en"]]
    config2 = [[[1, 2], "e"], [[1, 2], "n"], [[2, 3], "en"]]
    for [grades1, t1] in config1
      for g1 in grades1
        for [grades2, t2] in config2
          for g2 in grades2
            for o in [0, 2]
              mv = @mv [[g1, t1, 0], [g2, t2, o]]
              continue unless mv
              combinations.push [[mv[0]], [mv[1]]]
    mv = @mv [[3, "en", 0], [1, "e", 1]]
    combinations.push [[mv[0]], [mv[1]]]
    combinations.push [@mv([[1, "any", 0], [1, "any", 0]]), @mv([[2, "any", 0], [3, "any", 0]])]
    strings = for a in combinations
      a = (@space.mv_to_string(@space.mv(b)) for b in a)
      console.log @to_json a

class custom_test_runner_class extends test_runner_class
  # this converts input/output from and to strings for a compact and easier to manage test specification.
  is_mv: (a) -> Array.isArray(a) && a.every((a) -> Array.isArray(a) && 3 == a.length)
  is_mv_config: (a) -> Array.isArray(a) && a.every((a) -> Array.isArray(a) && 2 == a.length)
  is_mv_to_string: (a) -> if @is_mv_config a then @space.mv_to_string(@space.mv a) else a
  constructor: (space) ->
    super()
    @space = space
    @report_compact_failure_strings = (inp, exp, out) ->
      [((@is_mv_to_string a for a in inp).join " "), @is_mv_to_string(exp), out]

c3 = new sph_ga [1, 1, 1], conformal: true

# new test data is generated when needed.
#gen = new test_data_generator_class c3
#gen.display_combination_tests()

test_runner = new custom_test_runner_class c3

call_with_mv_string = (f) ->
  (...a) -> c3.mv_to_string(f.apply(c3, (c3.mv_from_string b for b in a)))

test_runner.execute [
  [
    ["ip", call_with_mv_string(c3.ip)]
    ["3e4", "4e1_3_4"]
    "0"
    ["2", "3"]
    "0"
    ["4", "5e1"]
    "20e1"
    ["6e2", "7"]
    "0"
    ["8e3", "9e1"]
    "0"
    ["2e2", "3e2"]
    "6"
    ["4e3", "5e1_2"]
    "0"
    ["4e2", "5e2_3"]
    "20e3"
    ["7e2", "8e4"]
    "0"
    ["2e1", "3e4_5"]
    "0"
    ["5e3", "6e1_4"]
    "0"
    ["8e2", "9e1_2_4"]
    "-72e1_4"
    ["3e5", "4e1"]
    "0"
    ["6e5", "7e1_3"]
    "0"
    ["9e5", "2e4"]
    "18"
    ["3e5", "4e5"]
    "0"
    ["5e4", "6e4_5"]
    "-30e4"
    ["8e4", "9e2_4"]
    "0"
    ["6e4_5", "7e2"]
    "0"
    ["9e4_5", "2e2_3"]
    "0"
    ["4e4_5", "5e4"]
    "0"
    ["7e4_5", "8e4_5"]
    "0"
    ["3e4_5", "4e3_4"]
    "0"
    ["6e4_5", "7e2_3_4"]
    "0"
    ["9e1_5", "2e3"]
    "0"
    ["4e3_5", "5e1_2"]
    "0"
    ["7e2_4", "8e5"]
    "0"
    ["2e1_5", "3e4_5"]
    "0"
    ["5e3_5", "6e1_4"]
    "0"
    ["7e2_4", "8e2_4"]
    "0"
    ["9e3_4", "2e1_2_5"]
    "0"
    ["4e1_3_5", "5e1"]
    "0"
    ["6e1 + 7e4", "8e4_5 + 9e2_3_5"]
    "-63e2_3 + -56e4"
  ]
  [
    ["id_bit_indices", ((...inp) -> (c3.id_bit_indices i for i in inp))]
    [0...(2 ** c3.n)]
    [[],[0],[1],[0,1],[2],[0,2],[1,2],[0,1,2],[3],[0,3],[1,3],[0,1,3],[2,3],[0,2,3],[1,2,3],
     [0,1,2,3],[4],[0,4],[1,4],[0,1,4],[2,4],[0,2,4],[1,2,4],[0,1,2,4],[3,4],[0,3,4],[1,3,4],
     [0,1,3,4],[2,3,4],[0,2,3,4],[1,2,3,4],[0,1,2,3,4]]
  ]
  [
    ["reverse", call_with_mv_string(c3.reverse)]
    "e1"
    "e1"
    "e1_2"
    "-1e1_2"
    "e1_2_3"
    "-1e1_2_3"
    "1 + e1"
    "1 + e1"
    "1 + e1 + e1_2"
    "1 + e1 + -1e1_2"
    "1 + e1 + e2 + e1_2 + e1_3 + e1_2_3"
    "1 + e1 + e2 + -1e1_2 + -1e1_3 + -1e1_2_3"
  ]
  [
    ["involute", call_with_mv_string(c3.involute)]
    "e1"
    "-1e1"
    "e1_2"
    "e1_2"
    "e1_2_3"
    "-1e1_2_3"
  ]
  [
    ["conjugate", call_with_mv_string(c3.conjugate)]
    "e1"
    "-1e1"
    "e1_2"
    "-1e1_2"
    "e1_2_3"
    "e1_2_3"
  ]
  [
    ["ep", call_with_mv_string(c3.ep)]
    ["3e4", "4e1_3_4"]
    "0"
    ["2", "3"]
    "0"
    ["4", "5e1"]
    "20e1"
    ["6e2", "7"]
    "42e2"
    ["8e3", "9e1"]
    "0"
    ["2e2", "3e2"]
    "6"
    ["4e3", "5e1_2"]
    "0"
    ["4e2", "5e2_3"]
    "20e3"
    ["7e2", "8e4"]
    "0"
    ["2e1", "3e4_5"]
    "0"
    ["5e3", "6e1_4"]
    "0"
    ["8e2", "9e1_2_4"]
    "-72e1_4"
    ["3e5", "4e1"]
    "0"
    ["6e5", "7e1_3"]
    "0"
    ["9e5", "2e4"]
    "18"
    ["3e5", "4e5"]
    "0"
    ["5e4", "6e4_5"]
    "-30e4"
    ["8e4", "9e2_4"]
    "0"
    ["6e4_5", "7e2"]
    "0"
    ["9e4_5", "2e2_3"]
    "0"
    ["4e4_5", "5e4"]
    "0"
    ["7e4_5", "8e4_5"]
    "0"
    ["3e4_5", "4e3_4"]
    "0"
    ["6e4_5", "7e2_3_4"]
    "0"
    ["9e1_5", "2e3"]
    "0"
    ["4e3_5", "5e1_2"]
    "0"
    ["7e2_4", "8e5"]
    "0"
    ["2e1_5", "3e4_5"]
    "0"
    ["5e3_5", "6e1_4"]
    "0"
    ["7e2_4", "8e2_4"]
    "0"
    ["9e3_4", "2e1_2_5"]
    "0"
    ["4e1_3_5", "5e1"]
    "0"
    ["6e1 + 7e4", "8e4_5 + 9e2_3_5"]
    "-63e2_3 + -56e4"
  ]
]
