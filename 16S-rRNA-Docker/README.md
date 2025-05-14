# 16S rRNA Microbial Analysis Docker Setup

This README provides instructions on how to set up and run the 16S rRNA microbial analysis project using Docker. The project utilizes R for data analysis and visualization of microbial community data.

## Project Structure

The project is organized as follows:

```
16S-rRNA-Docker
├── docker
│   ├── Dockerfile
│   └── r-requirements.txt
├── scripts
│   └── install_packages.R
├── .dockerignore
├── docker-compose.yml
├── README.md
└── .env
```

## Prerequisites

Before you begin, ensure you have the following installed on your machine:

- Docker
- Docker Compose

## Building the Docker Image

1. Navigate to the project directory:

   ```bash
   cd path/to/16S-rRNA-Docker
   ```

2. Build the Docker image using the provided Dockerfile:

   ```bash
   docker build -t 16s-rna-analysis docker/
   ```

## Running the Docker Container

You can run the Docker container using Docker Compose, which simplifies the process of managing multi-container Docker applications.

1. Start the services defined in `docker-compose.yml`:

   ```bash
   docker-compose up
   ```

2. This command will build the image (if not already built) and start the container. You can access the R environment within the container for analysis.

## Installing R Packages

The required R packages are listed in `docker/r-requirements.txt`. The `scripts/install_packages.R` script will install these packages when the container is built.

To manually install packages, you can run the following command inside the running container:

```R
source("scripts/install_packages.R")
```

## Stopping the Docker Container

To stop the running container, press `CTRL + C` in the terminal where Docker Compose is running. You can also stop it using:

```bash
docker-compose down
```

## Additional Information

For any questions or issues, please refer to the documentation or contact the project maintainers.

Happy analyzing!