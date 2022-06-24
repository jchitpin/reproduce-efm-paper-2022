## DESCRIPTION -----------------------------------------------------------------
# This script does the following (for wildtype/disease sphingolipid dataset):
# (1) Enumerate and compute EFM weights by cycle-history Markov chain method.
# (2) Export EFM weights.
# (3) Export EFM matrix.
# (4) Compute the ∑log10(squared reconstruction error) of Markov solution.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/manuscript-01-computations/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using ElementaryFluxModes # https://github.com/jchitpin/ElementaryFluxModes.jl
using Tables, PrettyTables
# ------------------------------------------------------------------------------

## LOAD STOICHIOMETRY AND FLUX DATA --------------------------------------------
stoich_wt = CSV.read("../data/stoich_wt.csv", Tables.matrix, header=false)
stoich_ad = CSV.read("../data/stoich_ad.csv", Tables.matrix, header=false)
fluxes_wt = vec(CSV.read("../data/fluxes_wt.csv", Tables.matrix, header=false))
fluxes_ad = vec(CSV.read("../data/fluxes_ad.csv", Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## COMPUTING EFM PROBABILITIES BY MARKOV CHAIN ---------------------------------
res_wt = steady_state_efm_distribution(stoich_wt, fluxes_wt)
res_ad = steady_state_efm_distribution(stoich_ad, fluxes_ad)
# ------------------------------------------------------------------------------

## EXPORTING EFM WEIGHTS -------------------------------------------------------
CSV.write(#
  "../data/efm_weights_wt_markov.csv",
  Tables.table(res_wt.w),
  header=false
)
CSV.write(#
  "../data/efm_weights_ad_markov.csv",
  Tables.table(res_ad.w),
  header=false
)
# ------------------------------------------------------------------------------

## EXPORTING EFM MATRIX --------------------------------------------------------
# EFM matrix (rows are reactions; columns are EFMs)
A_wt = reshape_efm_vector(res_wt.e, stoich_wt)
A_ad = reshape_efm_vector(res_ad.e, stoich_ad)
@assert(A_wt == A_ad)
CSV.write(#
  "../data/efm_matrix_sphingo.csv",
  Tables.table(A_wt),
  header=false
)
# ------------------------------------------------------------------------------

## LOG10 SQUARED RECONSTRUCTION ERROR ------------------------------------------
# ∑(log₁₀((Ax-v)²)) = -592.85 (wildtype)
error_wt = sum(filter!(!isinf, log10.((A_wt * res_wt.w .- fluxes_wt).^2)))
# ∑(log₁₀((Ax-v)²)) = -630.55 (disease)
error_ad = sum(filter!(!isinf, log10.((A_ad * res_ad.w .- fluxes_ad).^2)))
# ------------------------------------------------------------------------------

