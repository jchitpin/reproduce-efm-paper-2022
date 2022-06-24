# Any zero EFMs that are logged become -Inf. An EFM with a raw weight of 10^-10
# and zero are converted to a very small value (-900) and used in the following
# section to compute the frequency of zero EFM weights.
function aggregate_efm_values(#
  mc::Vector{Float64},
  combined,
  type::Symbol,
  fname::String
)
  @assert(length(combined) == 5, "Expected results from 5 objective function.")

  # Order data by Markov and optionally transform and remove infinities
  idx_mc = sortperm(mc)
  t(x) = x
  if type ∈ [:efms_raw_log, :efms_prop_log]
    f(x) = log10.(x)
    g1(x) = replace!(y -> isinf(y) ? -900 : y, x)
    g2(x) = replace!(y -> y <= -10 ? -900 : y, x)
    t(x) = g2(g1(f(x)))[idx_mc]
  elseif type ∈ [:efms_raw, :efms_prop]
    t(x) = x[idx_mc]
  else
    @assert(false, "type symbol value incorrect.")
  end

  # Indices to access combined
  if type == :efms_raw
    idx = 1
  elseif type == :efms_prop
    idx = 2
  elseif type == :efms_raw_log
    idx = 3
  elseif type == :efms_prop_log
    idx = 4
  end

  # EFM values
  combined_vals = vcat(#
    [#
      mc[idx_mc];
      t.(combined[1][idx]);
      t.(combined[2][idx]);
      t.(combined[3][idx]);
      t.(combined[4][idx]);
      t.(combined[5][idx]);
    ]...
  )

  # The following is hardcoded
  idx_vals = repeat(1:55, 22)
  efms_idx = repeat(idx_mc, 22)
  methods = [#
    repeat(["Markov"], 55);
    repeat(["L2_norm"], 55*2);
    repeat(["qp_max_spa"], 55*2);
    repeat(["lp_max_spa"], 55*8);
    repeat(["lp_min_spa"], 55*8);
    repeat(["milp"], 55);
  ]
  solvers = [#
    repeat(["Markov"], 55);
    repeat(["COSMO"], 55);
    repeat(["OSQP"], 55);
    repeat(["COSMO"], 55);
    repeat(["OSQP"], 55);
    repeat(["GLPK"], 55);
    repeat(["Gurobi"], 55);
    repeat(["SCIP"], 55);
    repeat(["CDDLib"], 55);
    repeat(["ECOS"], 55);
    repeat(["OSQP"], 55);
    repeat(["ProxSDP"], 55);
    repeat(["Tulip"], 55);
    repeat(["GLPK"], 55);
    repeat(["Gurobi"], 55);
    repeat(["SCIP"], 55);
    repeat(["CDDLib"], 55);
    repeat(["ECOS"], 55);
    repeat(["OSQP"], 55);
    repeat(["ProxSDP"], 55);
    repeat(["Tulip"], 55);
    repeat(["Gurobi"], 55);
  ]

  # EXPORT
  CSV.write(#
    fname,
    Tables.table(#
      hcat(#
        idx_vals, combined_vals, methods, solvers, efms_idx
      )
    ),
    header = ["x", "y", "class", "solver", "efm"]
  )
end

function aggregate_efm_zero_freq(mc::Vector{Float64}, combined)

  idx_mc = sortperm(mc)
  g1(x) = replace!(y -> isinf(y) ? -900 : y, x)
  g2(x) = replace!(y -> y <= -10 ? -900 : y, x)
  s(x) = findall(==(-900), g2(g1(x))[idx_mc])
  function u(x::Vector{Float64})
    y = zeros(length(x))
    y[s(x)] .= 1
    return y
  end

  bars = hcat(#
    sum(u.(combined[1][3])),
    sum(u.(combined[2][3])),
    sum(u.(combined[3][3])),
    sum(u.(combined[4][3])),
    sum(u.(combined[5][3]))
  )
  return bars ./ 21
end






