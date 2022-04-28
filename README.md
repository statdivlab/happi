
<!-- README.md is generated from README.Rmd. Please edit that file -->

`happi`: <span style="font-family:Arial; font-size:2em;"> a **H**ierarchical
**Ap**proach to **P**angenomics **I**nference</span>

`happi` is an `R` package for modeling gene presence that accounts for differences in genome quality. 

## Installation

    if (!require("remotes", quietly = TRUE))
        install.packages("remotes") # check that remotes is installed
    remotes::install_github("statdivlab/happi", build_vignettes = TRUE) # install happi using remotes and build vignettes
    library(happi)

## Usage

The vignettes provide examples of how to use `happi` and all its main
functions. You can follow the vignettes by running the following code in
`R`:

    utils::browseVignettes(package = "happi")

## Citation

If you use `happi` please cite our work:

An open-access preprint is available [here](https://www.biorxiv.org/content/10.1101/2022.04.26.489591v1).

## Issues/Requests

If you have any issues using our software or further questions please
submit an issue [here](https://github.com/statdivlab/happi/issues).
