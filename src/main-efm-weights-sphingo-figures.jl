## DESCRIPTION -----------------------------------------------------------------
# This script does the following (for wildtype/disease sphingolipid dataset):
# (1)  Compute EFM weights using LP, MILP, and QP formulations described
#      in literature using various solvers.
# (2a) Compute the log10(squared reconstruction error) for the Markov solution
# (2b) Compute the log10(squared reconstruction error) for optimization methods
# (3)  Export log10(squared reconstruction error) for optimization methods
# (4)  Export EFM weights from each method/solver.
# ------------------------------------------------------------------------------
#
## USER PARAMETERS -------------------------------------------------------------
# Set working directory of this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2022/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using Tables, CSV
using JuMP
using GLPK, Gurobi, SCIP, COSMO, OSQP, CDDLib, ECOS, ProxSDP, Tulip
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## LOAD RELEVANT DATA ----------------------------------------------------------
fluxes_wt = vec(CSV.read("../data/fluxes-wt.csv", Tables.matrix, header=false))
fluxes_ad = vec(CSV.read("../data/fluxes-ad.csv", Tables.matrix, header=false))
efms = CSV.read("../data/efm-matrix-sphingo.csv", Tables.matrix, header=false)
w_raw_wt = vec(CSV.read("../data/efm-weights-wt-raw-markov.csv", Tables.matrix, header=false))
w_raw_ad = vec(CSV.read("../data/efm-weights-ad-raw-markov.csv", Tables.matrix, header=false))
w_log_wt = vec(CSV.read("../data/efm-weights-wt-log-markov.csv", Tables.matrix, header=false))
w_log_ad = vec(CSV.read("../data/efm-weights-ad-log-markov.csv", Tables.matrix, header=false))
p_raw_wt = vec(CSV.read("../data/efm-proport-wt-raw-markov.csv", Tables.matrix, header=false))
p_raw_ad = vec(CSV.read("../data/efm-proport-ad-raw-markov.csv", Tables.matrix, header=false))
p_log_wt = vec(CSV.read("../data/efm-proport-wt-log-markov.csv", Tables.matrix, header=false))
p_log_ad = vec(CSV.read("../data/efm-proport-ad-log-markov.csv", Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## SETTING UP MATHEMATICAL OPTIMIZATION VARIABLES ------------------------------
# Renaming variables to follow convention: Ax=v
A = copy(efms)
v1 = fluxes_wt
v2 = fluxes_ad

# Solvers
milp_solvers, qp_solvers, lp_solvers = load_solvers()
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING --------------------------------------------------------
# Reference: Schwartz and Kanehisa (2006) DOI:10.1186/1471-2105-7-186
sk_wt = optimization_sk(A, v1, qp_solvers)
sk_ad = optimization_sk(A, v2, qp_solvers)

# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ---------------
# Reference: Orman et al. (2011) DOI:10.1016/j.jtbi.2010.11.042
or_wt = optimization_or(A, v1, qp_solvers)
or_ad = optimization_or(A, v2, qp_solvers)
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ------------------
# Reference: Rugen et al. (2012) DOI:10.1016/j.ymben.2012.01.009
ru_wt = optimization_ru(A, v1, lp_solvers)
ru_ad = optimization_ru(A, v2, lp_solvers)
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING LONGEST PATHWAYS -------------------------------
# Reference: Ren et al. (2020) DOI:10.1016/j.algal.2019.101767
re_wt = optimization_re(A, v1, lp_solvers)
re_ad = optimization_re(A, v2, lp_solvers)
# ------------------------------------------------------------------------------

# MIXED INTEGER PROGRAMMING MAXIMIZING SHORTEST ACTIVE PATH---------------------
# Reference: Nookaew et al. (2007) DOI:10.1002/bit.21339
no_wt = optimization_no(A, v1, milp_solvers)
no_ad = optimization_no(A, v2, milp_solvers)
# ------------------------------------------------------------------------------

# EXPLORATORY PLOTS ------------------------------------------------------------
# Correlation scatterplot
using Plots
l = @layout [a b c; d e f]
p1 = plot(w_log_wt, w_log_ad, seriestype = :scatter, titlefontsize = 6, title = "Markov", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
p2 = plot(sk_wt.efms_raw_log, vcat(sk_ad.efms_raw_log), seriestype = :scatter, titlefontsize = 6, title = "Min L2", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
p3 = plot(or_wt.efms_raw_log, or_ad.efms_raw_log, seriestype = :scatter, titlefontsize = 6, title = "Max qSPA", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
p4 = plot(ru_wt.efms_raw_log, ru_ad.efms_raw_log, seriestype = :scatter, titlefontsize = 6, title = "Max lSPA", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
p5 = plot(re_wt.efms_raw_log, re_ad.efms_raw_log, seriestype = :scatter, titlefontsize = 6, title = "Min lSPA", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
p6 = plot(no_wt.efms_raw_log, no_ad.efms_raw_log, seriestype = :scatter, titlefontsize = 6, title = "Min milAP", legend = :none, xlim = (-10, 1), ylim = (-10, 1))
plot(p1, p2, p3, p4, p5, p6, layout = l)

# FC sorted by original wildtype weights
l = @layout [a b c; d e f]
id = sortperm(w_raw_wt)
f(x,y,i) = log2.((x ./ y))[i]
p1 = plot(1:55, f(w_raw_ad, w_raw_wt, id), seriestype = :scatter, titlefontsize = 6, title = "Markov", legend = :none)
p2 = plot(1:55, f.(sk_ad.efms_raw, sk_wt.efms_raw, Ref(id)), seriestype = :scatter, titlefontsize = 6, title = "Min L2", legend = :none)
p3 = plot(1:55, f.(or_ad.efms_raw, or_wt.efms_raw, Ref(id)), seriestype = :scatter, titlefontsize = 6, title = "Max qSPA", legend = :none)
p4 = plot(1:55, f.(ru_ad.efms_raw, ru_wt.efms_raw, Ref(id)), seriestype = :scatter, titlefontsize = 6, title = "Max lSPA", legend = :none)
p5 = plot(1:55, f.(re_ad.efms_raw, re_wt.efms_raw, Ref(id)), seriestype = :scatter, titlefontsize = 6, title = "Min lSPA", legend = :none)
p6 = plot(1:55, f.(no_ad.efms_raw, no_wt.efms_raw, Ref(id)), seriestype = :scatter, titlefontsize = 6, title = "Min milAP", legend = :none)
plot(p1, p2, p3, p4, p5, p6, layout = l)

plot(1:55, f(w_raw_ad, w_raw_wt, id), seriestype = :scatter, legend = :none, markercolor = :red)
plot!(1:55, f.(sk_ad.efms_raw, sk_wt.efms_raw, Ref(id)), seriestype = :scatter, markercolor = :blue)
plot!(1:55, f.(or_ad.efms_raw, or_wt.efms_raw, Ref(id)), seriestype = :scatter, markercolor = :green)
plot!(1:55, f.(ru_ad.efms_raw, ru_wt.efms_raw, Ref(id)), seriestype = :scatter, markercolor = :orange)
plot!(1:55, f.(re_ad.efms_raw, re_wt.efms_raw, Ref(id)), seriestype = :scatter, markercolor = :purple)
plot!(1:55, f.(no_ad.efms_raw, no_wt.efms_raw, Ref(id)), seriestype = :scatter, markercolor = :yellow)

## Favourite plot for now - Fold change ordered by sorted order of flux explained
# Only one solver for each method (one with the lowest error across wildtype/disease)
f(x,y,i) = log2.((x ./ y))[i]
g(x,y) = log10.(sort(x .* y))/sum(x .* y)

p1 = plot(g(w_raw_wt, vec(sum(A, dims=1))), f(w_raw_ad, w_raw_wt, id), seriestype = :scatter, legend = :none, markercolor = :red, xlim = (-7, 0), ylim = (-20, 10), xlabel = "EFMs sorted by how much flux they explain log₁₀(%)", ylabel = "Fold change log₂(AD/WT)")
plot!(g.(sk_wt.efms_raw[[1]], Ref(vec(sum(A, dims=1)))), f.(sk_ad.efms_raw[[1]], sk_wt.efms_raw[[1]], sortperm.(sk_wt.efms_raw[[1]])), seriestype = :scatter, markercolor = :blue)
plot!(g.(or_wt.efms_raw[[1]], Ref(vec(sum(A, dims=1)))), f.(or_ad.efms_raw[[1]], or_wt.efms_raw[[1]], sortperm.(or_wt.efms_raw[[1]])), seriestype = :scatter, markercolor = :green)
plot!(g.(ru_wt.efms_raw[[3]], Ref(vec(sum(A, dims=1)))), f.(ru_ad.efms_raw[[3]], ru_wt.efms_raw[[3]], sortperm.(ru_wt.efms_raw[[3]])), seriestype = :scatter, markercolor = :orange)
plot!(g.(re_wt.efms_raw[[3]], Ref(vec(sum(A, dims=1)))), f.(re_ad.efms_raw[[3]], re_wt.efms_raw[[3]], sortperm.(re_wt.efms_raw[[3]])), seriestype = :scatter, markercolor = :purple)
plot!(g.(no_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(no_ad.efms_raw, no_wt.efms_raw, sortperm.(no_wt.efms_raw)), seriestype = :scatter, markercolor = :yellow)
plot!([-7; 0], [2.5; 2.5], lw=1, lc=:black, legend=false, line=(:dash))
plot!([-7; 0], [-2.5; -2.5], lw=1, lc=:black, legend=false, line=(:dash))

p2 = plot(g(w_raw_wt, vec(sum(A, dims=1))), f(w_raw_ad, w_raw_wt, id), seriestype = :scatter, legend = :none, markercolor = :red, xlim = (-7, 0), ylim = (-20, 10), xlabel = "EFMs sorted by how much flux they explain log₁₀(%)", ylabel = "Fold change log₂(AD/WT)")
plot!(g.(sk_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(sk_ad.efms_raw, sk_wt.efms_raw, sortperm.(sk_wt.efms_raw)), seriestype = :scatter, markercolor = :blue)
plot!(g.(or_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(or_ad.efms_raw, or_wt.efms_raw, sortperm.(or_wt.efms_raw)), seriestype = :scatter, markercolor = :green)
plot!(g.(ru_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(ru_ad.efms_raw, ru_wt.efms_raw, sortperm.(ru_wt.efms_raw)), seriestype = :scatter, markercolor = :orange)
plot!(g.(re_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(re_ad.efms_raw, re_wt.efms_raw, sortperm.(re_wt.efms_raw)), seriestype = :scatter, markercolor = :purple)
plot!(g.(no_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(no_ad.efms_raw, no_wt.efms_raw, sortperm.(no_wt.efms_raw)), seriestype = :scatter, markercolor = :yellow)
plot!([-9; 0.1], [2.5; 2.5], lw=1, lc=:black, legend=false, line=(:dash))
plot!([-9; 0.1], [-2.5; -2.5], lw=1, lc=:black, legend=false, line=(:dash))


# ------------------------------------------------------------------------------




## CORRELATING WILDTYPE AND DISEASE EFM WEIGHTS --------------------------------
correlate_efm_values(#
  w_log_wt,
  w_log_ad,
  combined_wt,
  combined_ad,
  :efms_prop_log,
  "../data/scatterplot-correlate-efm-weights-log.csv"
)
# ------------------------------------------------------------------------------

