# Contents of /16S-rRNA-Docker/16S-rRNA-Docker/scripts/install_packages.R

# This script installs the R packages listed in r-requirements.txt

# Read the package names from the requirements file
packages <- readLines("r-requirements.txt")

# Install the packages
install.packages(packages, repos = "http://cran.r-project.org")