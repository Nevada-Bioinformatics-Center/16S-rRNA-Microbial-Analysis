FROM rocker/rstudio:4.4.2

# Set the working directory
WORKDIR /usr/local/src/microbiome-analysis

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the R requirements file
COPY docker/r-requirements.txt .

# Install R packages
COPY scripts/install_packages.R .
RUN Rscript install_packages.R

# Clean up
RUN rm install_packages.R
