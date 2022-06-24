# Benchmarking EFM weights by cycle-history Markov chain versus optimization-based methods

This repository contains the scripts necessary to regenerate all figures and
results from the manuscript by J. G. Chitpin and T. J. Perkins.

## Hardware requirements

Any architecture that can run the Julia programming language.


## Software requirements

The following is required:

* Julia (minimum version 1.6).
* A Bash terminal to run the scripts/workflows.
* TeX distribution (like TeX Live or MiKTeX) to compile figures.

## Installation

Download the repository and install all necessary Julia packages by running the
following commands in your desired installation directory.

`$ cd /home/<username>/<directory>/`  
`$ git clone jchitpin/reproduce-efm-paper-2022`  
`$ cd reproduce-efm-paper-2022/src/`  
`$ julia install_julia_packages.jl` # or run line by line in Julia REPL

## Workflow to reproduce results

The scripts in the following subsections should be run in order to regenerate
the intermediate data files.

### Figure 1 and 2

1. ``



### Figure 3


### Figure 4

1. `$ julia main_sphingo_network_validation.jl`



## Reference



## Acknowledgements

This work was supported in part by [grant] from [agency]. J.G.C. was
supported by an NSERC CREATE Matrix Metabolomics Scholarship and
an NSERC Alexander Graham Bell Canada Graduate Scholarship.




