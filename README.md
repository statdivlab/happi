
<!-- README.md is generated from README.Rmd. Please edit that file -->

<span style="font-family:Arial; font-size:2em;"> a **H**ierarchical
**Ap**proach to **P**angenomics **I**nference</span>

`happi` is an `R` package for modeling gene presence.

## Installation

    if (!require("devtools", quietly = TRUE))
        install.packages("devtools") # check that devtools is installed
    devtools::install_github("statdivlab/happi", build_vignettes = TRUE) # install happi using devtools
    library(happi)

## Usage

The vignettes provide examples of how to use `happi` and all its main
functions through the `R` interactive session. You can follow the
vignettes by running the following code in `R`:

    utils::browseVignettes(package = "happi")

An example snakemake workflow of `happi`’s usage has been made available
under the `workflows/` folder of this github directory. To run the
example workflow you’ll need to install snakemake:

    conda install -n base -c conda-forge mamba
    conda activate base
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

    conda activate snakemake

Once your snakemake has been installed in a conda environment and
activated you can run the following to execute the workflow:

    snakemake -s Snakefile --cores 6

## Citation

If you use `happi` please cite our work:

An open-access preprint is available here.

## Issues/Requests

If you have any issues using our software or further questions please
submit an issue [here](https://github.com/statdivlab/happi/issues).
