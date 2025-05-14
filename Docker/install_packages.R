#!/usr/bin/env Rscript

# R script to install requirements for exercises -------------------------------

## a vector of packages to install (edit in this section) ----------------------
### packages could be either on CRAN or bioconductor

# force compilation from source for tidytree
# binary installation is only available for earlier versions



pkgs <- c("devtools", "remotes",
  "dplyr",
  "dada2", "ggplot2",
  "phyloseq", "vegan",
  "lubridate", "dplyr",
  "tidyverse", "ggpubr",
  "rstatix", "stringr"
)

## install Bioconductor --------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## install and check package loading -------------------------------------------
for (pkg in basename(pkgs)) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)

  if (!library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(
      paste0(
        "Installation of package ",
        pkg,
        " exited with non-zero exit status"
      ),
      stdout()
    )
    quit(status = 1, save = "no")
  }
}


# From repo
devtools::install_github("thomasp85/patchwork")
remotes::install_github("gavinsimpson/ggvegan")

