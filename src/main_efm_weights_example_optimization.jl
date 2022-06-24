## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Compute EFM weights in the two example networks by optimization-based
#     methods across different solvers.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/manuscript-01-computations/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using CSV, Tables
using JuMP # interface for mathematical optimization (solvers are listed below)
using GLPK, Gurobi, SCIP, COSMO, OSQP, CDDLib, ECOS, ProxSDP, Tulip
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## SETTING UP MATHEMATICAL OPTIMIZATION VARIABLES ------------------------------
# Goal is to create system of linear EFM weight equations for example networks
# 1 and 2 with the following format: A * w == v

# Enumerate EFMs using ElementaryFluxModes and reshape into EFM matrix
S1 = [#
    -1 -1  0  0  0  0  0  0  1
     1  0 -1  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0
     0  0  1  1 -1 -1  0  0  0
     0  0  0  0  1  0 -1  0  0
     0  0  0  0  0  1  0 -1  0
     0  0  0  0  0  0  1  1 -1
]

# Steady state fluxes
v1 = [2, 2, 2, 2, 2, 2, 2, 2, 4]

# Enumerate EFMs, compute their probabilities and weights
res1 = steady_state_efm_distribution(S1, v1)

# EFM matrix (rows are reactions; columns are EFMs)
A1 = reshape_efm_vector(res1.e, S1)

# Stoichiometry matrix (rows are metabolites; cols are reactions)
S2 = [#
    -1  0  0  0  0  0  0  0  0  0  1
     1 -1  1 -1  0  0  0  0  0  0  0
     0  1 -1  0 -1  1  0  0  0  0  0
     0  0  0  1  0  0 -1  0  0  0  0
     0  0  0  0  1 -1  1 -1  1 -1  0
     0  0  0  0  0  0  0  0  0  1 -1
     0  0  0  0  0  0  0  1 -1  0  0
]

# Steady state fluxes
v2 = [3, 2, 1, 2, 3, 2, 2, 1, 1, 3, 3]

# Enumerate EFMs, compute their probabilities and weights
res2 = steady_state_efm_distribution(S2, v2)

# EFM matrix (rows are reactions; columns are EFMs)
A2 = reshape_efm_vector(res2.e, S2)

# Solvers
milp_solvers, qp_solvers, lp_solvers = load_solvers()
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING --------------------------------------------------------
# Reference: Schwartz and Kanehisa (2006) DOI:10.1186/1471-2105-7-186
efms_1_sk = Vector{Vector{Float64}}()
efms_2_sk = Vector{Vector{Float64}}()
for solver in qp_solvers
  push!(efms_1_sk, jump_qp_l2norm(A1, v1, solver))
  push!(efms_2_sk, jump_qp_l2norm(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ---------------
# Reference: Orman et al. (2011) DOI:10.1016/j.jtbi.2010.11.042
efms_1_or = Vector{Vector{Float64}}()
efms_2_or = Vector{Vector{Float64}}()
for solver in qp_solvers
  push!(efms_1_or, jump_qp_max_shortest_paths(A1, v1, solver))
  push!(efms_2_or, jump_qp_max_shortest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ------------------
# Reference: Rugen et al. (2012) DOI:10.1016/j.ymben.2012.01.009
efms_1_ru = Vector{Vector{Float64}}()
efms_2_ru = Vector{Vector{Float64}}()
for solver in lp_solvers
  push!(efms_1_ru, jump_lp_max_shortest_paths(A1, v1, solver))
  push!(efms_2_ru, jump_lp_max_shortest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING LONGEST PATHWAYS -------------------------------
# Reference: Ren et al. (2020) DOI:10.1016/j.algal.2019.101767
efms_1_re = Vector{Vector{Float64}}()
efms_2_re = Vector{Vector{Float64}}()
for solver in lp_solvers
  push!(efms_1_re, jump_lp_max_longest_paths(A1, v1, solver))
  push!(efms_2_re, jump_lp_max_longest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# MIXED INTEGER PROGRAMMING MAXIMIZING SHORTEST ACTIVE PATH---------------------
# Reference: Nookaew et al. (2007) DOI:10.1002/bit.21339
efms_1_no = Vector{Vector{Float64}}()
efms_2_no = Vector{Vector{Float64}}()
for solver in milp_solvers
  push!(efms_1_no, jump_milp_max_active_paths(A1, v1, solver))
  push!(efms_2_no, jump_milp_max_active_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

