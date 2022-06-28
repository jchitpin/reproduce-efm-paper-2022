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

# EXPORT LOG10 RECONSTRUCTION ERRORS -------------------------------------------
# Markov precision over individual fluxes
error_wt_ad_mc_ind = vec(# error is -6.80 and -7.01 for wildtype/disease
  CSV.read(#
    "../data/efm-error-ind-wt-ad-markov.csv", Tables.matrix, header=false
  )
)

# Markov precision over total fluxes
error_wt_ad_mc_tot = vec(# error is -15.65 and -15.35 for wildtype/disease
  CSV.read(#
    "../data/efm-error-total-wt-ad-markov.csv", Tables.matrix, header=false
  )
)

# Optimization-based error (wildtype)
combined_error_wt_ind = hcat(#
  [#
    [-sk_wt.efms_raw_error_ind; repeat([""], 7)],
    [-or_wt.efms_raw_error_ind; repeat([""], 7)],
    [""; -ru_wt.efms_raw_error_ind],
    [""; -re_wt.efms_raw_error_ind],
    [""; ""; -no_wt.efms_raw_error_ind; repeat([""], 6)]
  ]...
)
combined_error_wt_ind = permutedims(string.(combined_error_wt_ind))
combined_error_wt_ind = hcat(string.(collect(1:5)), combined_error_wt_ind)
CSV.write(#
  "../data/error-wt-ind.csv", Tables.table(combined_error_wt_ind), header=false
)

combined_error_wt_tot = hcat(#
  [#
    [-sk_wt.efms_raw_error_tot; repeat([""], 7)],
    [-or_wt.efms_raw_error_tot; repeat([""], 7)],
    [""; -ru_wt.efms_raw_error_tot],
    [""; -re_wt.efms_raw_error_tot],
    [""; ""; -no_wt.efms_raw_error_tot; repeat([""], 6)]
  ]...
)
combined_error_wt_tot = permutedims(string.(combined_error_wt_tot))
combined_error_wt_tot = hcat(string.(collect(1:5)), combined_error_wt_tot)
CSV.write(#
  "../data/error-wt-tot.csv", Tables.table(combined_error_wt_tot), header=false
)

# Optimization-based error (disease)
combined_error_ad_ind = hcat(#
  [#
    [-sk_ad.efms_raw_error_ind; repeat([""], 7)],
    [-or_ad.efms_raw_error_ind; repeat([""], 7)],
    [""; -ru_ad.efms_raw_error_ind],
    [""; -re_ad.efms_raw_error_ind],
    [""; ""; -no_ad.efms_raw_error_ind; repeat([""], 6)]
  ]...
)
combined_error_ad_ind = permutedims(string.(combined_error_ad_ind))
combined_error_ad_ind = hcat(string.(collect(1:5)), combined_error_ad_ind)
CSV.write(#
  "../data/error-ad-ind.csv", Tables.table(combined_error_ad_ind), header=false
)
combined_error_ad_tot = hcat(#
  [#
    [-sk_ad.efms_raw_error_tot; repeat([""], 7)],
    [-or_ad.efms_raw_error_tot; repeat([""], 7)],
    [""; -ru_ad.efms_raw_error_tot],
    [""; -re_ad.efms_raw_error_tot],
    [""; ""; -no_ad.efms_raw_error_tot; repeat([""], 6)]
  ]...
)
combined_error_ad_tot = permutedims(string.(combined_error_ad_tot))
combined_error_ad_tot = hcat(string.(collect(1:5)), combined_error_ad_tot)
CSV.write(#
  "../data/error-ad-tot.csv", Tables.table(combined_error_ad_tot), header=false
)
# ------------------------------------------------------------------------------

## EXPORT LOG10 EFM WEIGHTS ACROSS ALL METHODS ---------------------------------
combined_wt = [sk_wt, or_wt, ru_wt, re_wt, no_wt]
combined_ad = [sk_ad, or_ad, ru_ad, re_ad, no_ad]
aggregate_efm_values(#
  w_log_wt,
  combined_wt,
  :efms_raw_log,
  "../data/scatterplot-total-efm-weights-wt-log.csv"
)
aggregate_efm_values(#
  w_log_ad,
  combined_ad,
  :efms_raw_log,
  "../data/scatterplot-total-efm-weights-ad-log.csv"
)
# ------------------------------------------------------------------------------

## EXPORT LOG10 EFM PROPORTIONS ACROSS ALL METHODS -----------------------------
combined_wt = [sk_wt, or_wt, ru_wt, re_wt, no_wt]
combined_ad = [sk_ad, or_ad, ru_ad, re_ad, no_ad]
aggregate_efm_values(#
  p_log_wt,
  combined_wt,
  :efms_prop_log,
  "../data/scatterplot-total-efm-proport-wt-log.csv"
)
aggregate_efm_values(#
  p_log_ad,
  combined_ad,
  :efms_prop_log,
  "../data/scatterplot-total-efm-proport-ad-log.csv"
)
# ------------------------------------------------------------------------------

## EXPORT FREQUENCY OF ZERO EFM WEIGHTS ----------------------------------------
bars_wt = aggregate_efm_zero_freq(w_log_wt, combined_wt)
CSV.write(#
  "../data/barplot-zero-weight-freq-wt.csv",
  Tables.table(hcat(1:size(bars_wt)[1], bars_wt)),
  header = ["Label", "S1", "S2", "S3", "S4", "S5"],
  delim = '\t'
)

bars_ad = aggregate_efm_zero_freq(w_log_ad, combined_ad)
CSV.write(#
  "../data/barplot-zero-weight-freq-ad.csv",
  Tables.table(hcat(1:size(bars_ad)[1], bars_ad)),
  header = ["Label", "S1", "S2", "S3", "S4", "S5"],
  delim = '\t'
)
# ------------------------------------------------------------------------------

## ZERO WEIGHT PROBABILITY -----------------------------------------------------
# Mean probability of all optimization methods/solvers assigning a zero weight
mean_zeros_wt = sum(sum(bars_wt, dims=2)) / size(bars_wt,1) # 0.40
mean_zeros_ad = sum(sum(bars_ad, dims=2)) / size(bars_ad,1) # 0.44
# ------------------------------------------------------------------------------

## PERCENTAGE OF FLUX EXPLAINED BY EACH EFM WEIGHT -----------------------------
# First element is Markov solution, rest are the other methods/solvers
num_efms_explain_flux_wt = efm_flux_percentage(w_raw_wt, combined_wt, A, 0.95)
num_efms_explain_flux_ad = efm_flux_percentage(w_raw_ad, combined_ad, A, 0.95)
# Mean of 13.333 EFMs explain 95% of fluxes
sum(num_efms_explain_flux_wt[2:end]) / length(num_efms_explain_flux_wt[2:end])
# Mean of 13 EFMs explain 95% of fluxes
sum(num_efms_explain_flux_ad[2:end]) / length(num_efms_explain_flux_ad[2:end])
# ------------------------------------------------------------------------------

## DIFFERENCE BETWEEN MARKOV WILDTYPE AND DISEASE EFM WEIGHTS ------------------
# Run main-efm-weights-sphingo-markov.jl first to get this text file
dat = CSV.read("../data/table-markov-sphingo.csv", Tables.matrix, header=true)

mc_top_wt = efm_flux_percentage_top(w_raw_wt, A, 0.95)
mc_top_ad = efm_flux_percentage_top(w_raw_wt, A, 0.95)
mc_top_wt == mc_top_ad

w, x, y, z = groupsort(dat, 4, 6)
sk = log2.(sk_ad.efms_raw[1]./sk_wt.efms_raw[1])
or = log2.(or_ad.efms_raw[1]./or_wt.efms_raw[1])
ru = log2.(ru_ad.efms_raw[3]./ru_wt.efms_raw[3])
re = log2.(re_ad.efms_raw[3]./re_wt.efms_raw[3])
no = log2.(no_ad.efms_raw[1]./no_wt.efms_raw[1])

using Plots
plot(x, [y, sk[w], or[w], ru[w], re[w], no[w]], seriestype = :scatter, title = "log₂(ad/wt)")

findall(in(mc_top_wt), w)




plot(dat[:,1], dat[:,4][idx1], seriestype = :scatter, title = "log₂(ad/wt)")
plot(dat[:,1], dat[:,5][idx2], seriestype = :scatter, title = "ad-wt")



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






