# Uniquely identifying EFM weights by cycle-history Markov chain

This repository contains the scripts necessary to regenerate all figures and
results from the manuscript by J. G. Chitpin and T. J. Perkins. Certain values,
such as mean reconstruction error for the Markovian weights, are not exported
as a text file but can be viewed by interactively running the scripts. The
beginning of each script will list what is computed/exported.

## Requirements

The following is required:

* Julia (minimum version 1.6).
* Gurobi (must be version 9.12; free academic license available)
* A shell to run the scripts/workflows.
* TeX distribution (like TeX Live or MiKTeX) to compile figures.


## Installation

Download the repository and install all necessary Julia packages by running the
following commands in your desired installation directory.

1. `$ cd /home/<username>/<directory>/`  
2. `$ git clone jchitpin/reproduce-efm-paper-2022`  
3. `$ cd reproduce-efm-paper-2022/src/`  
4. `$ julia install-julia-packages.jl # or run line by line in Julia REPL`


## Workflow to reproduce results

The scripts in the following subsections should be run in order to regenerate
the intermediate data files. All scripts should be run in their current working
directory (`reproduce-efm-paper-2022/src/`).

### Figure 1 and 2

1. `$ julia main-efm-weights-example-markov.jl`
2. `$ julia main-efm-weights-example-optimization.jl`

Figures are generated via:

1. `$ sh figure-01.sh`
2. `$ sh figure-02.sh`

### Figure 3

Figure is generated via:

1. `sh figure-03.sh`

### Figure 4

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. FILL
2. FILL

### Supplementary Figure 1

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. `sh supplementary-01.sh`

### Supplementary Figure 2

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. `sh supplementary-02.sh`

### Supplementary Table 1

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`

Table is generated via:

1. `sh supplementary-03.sh`

### Supplementary Figure 4

1. `$ matlab -nodisplay -nosplash -nodesktop -r benchmark_efm_matlab.m
2. `$ julia benchmark-efm-julia.jl

Figure is generated via:

1. `sh supplementary-04.sh`

## Reference

## Acknowledgements

This work was supported in part by [grant] from [agency]. J.G.C. was
supported by an NSERC CREATE Matrix Metabolomics Scholarship and
an NSERC Alexander Graham Bell Canada Graduate Scholarship.




