## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Load stoichiometry matrices generated/benchmarked in MATLAB with EFMs
#     computed using FluxModeCalculator.
# (2) Benchmark the same stoichiometry matrices with the Julia code.
# (3) Export the data for plotting.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2022/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using ElementaryFluxModes
using Tables, CSV
using BenchmarkTools
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## LOAD SIMULATED STOICHIOMETRY MATRICES AND FLUXES ----------------------------
S01 = CSV.read("../data/stoich-20-30.csv", Tables.matrix, header=false)
S02 = CSV.read("../data/stoich-20-35.csv", Tables.matrix, header=false)
S03 = CSV.read("../data/stoich-20-40.csv", Tables.matrix, header=false)
S04 = CSV.read("../data/stoich-20-45.csv", Tables.matrix, header=false)
S05 = CSV.read("../data/stoich-20-50.csv", Tables.matrix, header=false)
S06 = CSV.read("../data/stoich-20-55.csv", Tables.matrix, header=false)
S07 = CSV.read("../data/stoich-20-60.csv", Tables.matrix, header=false)
S08 = CSV.read("../data/stoich-20-65.csv", Tables.matrix, header=false)
S09 = CSV.read("../data/stoich-20-70.csv", Tables.matrix, header=false)
S10 = CSV.read("../data/stoich-20-75.csv", Tables.matrix, header=false)
S11 = CSV.read("../data/stoich-20-80.csv", Tables.matrix, header=false)
v01 = vec(CSV.read("../data/ss-flux-20-30.csv", Tables.matrix, header=false))
v02 = vec(CSV.read("../data/ss-flux-20-35.csv", Tables.matrix, header=false))
v03 = vec(CSV.read("../data/ss-flux-20-40.csv", Tables.matrix, header=false))
v04 = vec(CSV.read("../data/ss-flux-20-45.csv", Tables.matrix, header=false))
v05 = vec(CSV.read("../data/ss-flux-20-50.csv", Tables.matrix, header=false))
v06 = vec(CSV.read("../data/ss-flux-20-55.csv", Tables.matrix, header=false))
v07 = vec(CSV.read("../data/ss-flux-20-60.csv", Tables.matrix, header=false))
v08 = vec(CSV.read("../data/ss-flux-20-65.csv", Tables.matrix, header=false))
v09 = vec(CSV.read("../data/ss-flux-20-70.csv", Tables.matrix, header=false))
v10 = vec(CSV.read("../data/ss-flux-20-75.csv", Tables.matrix, header=false))
v11 = vec(CSV.read("../data/ss-flux-20-80.csv", Tables.matrix, header=false))

E01 = CSV.read("../data/stoich-20-30-efm-matrix.csv", Tables.matrix, header=false)
E02 = CSV.read("../data/stoich-20-35-efm-matrix.csv", Tables.matrix, header=false)
E03 = CSV.read("../data/stoich-20-40-efm-matrix.csv", Tables.matrix, header=false)
E04 = CSV.read("../data/stoich-20-45-efm-matrix.csv", Tables.matrix, header=false)
E05 = CSV.read("../data/stoich-20-50-efm-matrix.csv", Tables.matrix, header=false)
E06 = CSV.read("../data/stoich-20-55-efm-matrix.csv", Tables.matrix, header=false)
# ------------------------------------------------------------------------------

## BENCHMARKING JULIA CODE TO ENUMERATE EFMS AND ASSIGN WEIGHTS ----------------
# Set benchmarking parameters for 100 samples
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 12_000
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
@benchmark steady_state_efm_distribution(S01, v01)
@benchmark steady_state_efm_distribution(S03, v03)
@benchmark steady_state_efm_distribution(S05, v05)
@benchmark steady_state_efm_distribution(S06, v06)  # 360 seconds
#@benchmark steady_state_efm_distribution(S07, v07) # infeasible
#@benchmark steady_state_efm_distribution(S08, v08) # infeasible
#@benchmark steady_state_efm_distribution(S09, v09) # infeasible
#@benchmark steady_state_efm_distribution(S10, v10) # infeasible
#@benchmark steady_state_efm_distribution(S11, v11) # infeasible
# ------------------------------------------------------------------------------

## CONFIRMING JULIA CODE ENUMERATES SAME NUMBER OF EFMS ------------------------
R01 = steady_state_efm_distribution(S01, v01)
R02 = steady_state_efm_distribution(S02, v02)
R03 = steady_state_efm_distribution(S03, v03)
R04 = steady_state_efm_distribution(S04, v04)
R05 = steady_state_efm_distribution(S05, v05)
R06 = steady_state_efm_distribution(S06, v06)
#R07 = steady_state_efm_distribution(S07, v07)
#R08 = steady_state_efm_distribution(S08, v08)
#R09 = steady_state_efm_distribution(S09, v09)
#R10 = steady_state_efm_distribution(S10, v10)
#R11 = steady_state_efm_distribution(S11, v11)

# It so happens that FluxModeCalculator returns duplicate EFMs while
# the Julia code returns the exact number of EFMs. However, the Julia code
# fails to scale beyond approximately 2000 EFMs. The memory requirement to
# generate the cycle-history Markov chain increases exponentially with EFM
# count and the code seems to fail when computing the steady state
# probabilities using the QuantEcon package
idx_intersect_01 = compare_efms(E01, R01.e, S01);
idx_intersect_02 = compare_efms(E02, R02.e, S02);
idx_intersect_03 = compare_efms(E03, R03.e, S03);
idx_intersect_04 = compare_efms(E04, R04.e, S04);
idx_intersect_05 = compare_efms(E05, R05.e, S05);
idx_intersect_06 = compare_efms(E06, R06.e, S06);
#idx_intersect_07 = compare_efms(E07, R07.e, S07);
#idx_intersect_08 = compare_efms(E08, R08.e, S08);
#idx_intersect_09 = compare_efms(E09, R09.e, S09);
#idx_intersect_10 = compare_efms(E10, R10.e, S10);
#idx_intersect_11 = compare_efms(E11, R11.e, S11);
# ------------------------------------------------------------------------------

