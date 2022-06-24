## DESCRIPTION -----------------------------------------------------------------
# This script does the following (for wildtype/disease sphingolipid dataset):
# (1) Enumerate and compute EFM weights by cycle-history Markov chain method.
# (2) Export EFM weights.
# (3) Export EFM matrix.
# (4) Compute and export log10(∑(squared reconstruction error)/|v|)
#     for the Markov solution.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2022/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using ElementaryFluxModes
using Tables, CSV
# ------------------------------------------------------------------------------

## LOAD STOICHIOMETRY AND FLUX DATA --------------------------------------------
stoich_wt = CSV.read("../data/stoich-corrected.csv", Tables.matrix, header=false)
stoich_ad = CSV.read("../data/stoich-corrected.csv", Tables.matrix, header=false)
fluxes_wt = vec(CSV.read("../data/fluxes-wt.csv", Tables.matrix, header=false))
fluxes_ad = vec(CSV.read("../data/fluxes-ad.csv", Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## COMPUTING EFM PROBABILITIES BY MARKOV CHAIN ---------------------------------
res_wt = steady_state_efm_distribution(stoich_wt, fluxes_wt)
res_ad = steady_state_efm_distribution(stoich_ad, fluxes_ad)
# ------------------------------------------------------------------------------

## EXPORTING EFM PROPORTIONS ---------------------------------------------------
# Raw weights
CSV.write(#
  "../data/efm-proport-wt-raw-markov.csv",
  Tables.table(res_wt.p),
  header=false
)
CSV.write(#
  "../data/efm-proport-ad-raw-markov.csv",
  Tables.table(res_ad.p),
  header=false
)
# Log10(weights)
CSV.write(#
  "../data/efm-proport-wt-log-markov.csv",
  Tables.table(log10.(res_wt.p)),
  header=false
)
CSV.write(#
  "../data/efm-proport-ad-log-markov.csv",
  Tables.table(log10.(res_ad.p)),
  header=false
)
# ------------------------------------------------------------------------------

## EXPORTING EFM WEIGHTS -------------------------------------------------------
# Raw weights
CSV.write(#
  "../data/efm-weights-wt-raw-markov.csv",
  Tables.table(res_wt.w),
  header=false
)
CSV.write(#
  "../data/efm-weights-ad-raw-markov.csv",
  Tables.table(res_ad.w),
  header=false
)
# Log10(weights)
CSV.write(#
  "../data/efm-weights-wt-log-markov.csv",
  Tables.table(log10.(res_wt.w)),
  header=false
)
CSV.write(#
  "../data/efm-weights-ad-log-markov.csv",
  Tables.table(log10.(res_ad.w)),
  header=false
)
# ------------------------------------------------------------------------------

## EXPORTING EFM MATRIX --------------------------------------------------------
# EFM matrix (rows are reactions; columns are EFMs)
A_wt = reshape_efm_vector(res_wt.e, stoich_wt)
A_ad = reshape_efm_vector(res_ad.e, stoich_ad)
@assert(A_wt == A_ad)
CSV.write(#
  "../data/efm-matrix-sphingo.csv",
  Tables.table(A_wt),
  header=false
)
# ------------------------------------------------------------------------------

## LOG10 SQUARED RECONSTRUCTION ERROR ------------------------------------------
# log₁₀(∑((Ax-v)²))/|v| = -6.80 (wildtype)
error_wt = log10(sum(((A_wt * res_wt.w .- fluxes_wt).^2)) / length(fluxes_wt))
# log₁₀(∑((Ax-v)²))/|v| = -7.01 (disease)
error_ad = log10(sum(((A_ad * res_ad.w .- fluxes_ad).^2)) / length(fluxes_ad))
CSV.write(#
  "../data/efm-error-wt-ad-markov.csv",
  Tables.table([error_wt, error_ad]),
  header=false
)
# ------------------------------------------------------------------------------

