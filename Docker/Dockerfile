FROM continuumio/miniconda3
ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

ENV NB_USER=$NB_USER \
    NB_UID=$NB_UID \
    NB_GID=$NB_GID \
    HOME=/home/$NB_USER


USER root
RUN useradd --non-unique --create-home --shell /bin/bash --gid "${NB_GID}" --uid ${NB_UID} "sagemaker-user"
# R pre-requisites

RUN apt-get update --yes && \
    apt-get install --yes software-properties-common lsb-release && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get install --yes --no-install-recommends \
    cmake \
    fonts-dejavu \
    unixodbc \
    unixodbc-dev \
    r-cran-rodbc \
    gfortran \
    awscli \
    gcc \
    r-base \
    r-base-dev \
    build-essential \ 
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev  \
    libjpeg-dev \
    libmagick++-dev \
    sudo \
    libcairo2-dev -y && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install jupyterlab && python3 -m pip install --upgrade pip &&  python3 -m pip install --upgrade urllib3==1.26.6
    
# R packages including IRKernel which gets installed globally.
# r-e1071: dependency of the caret R package
RUN conda config --add channels conda-forge &&\
    conda update --all --override-channels -c conda-forge --yes && \
    conda install mamba && \
    mamba install --quiet --yes \
    'boto3' \
    # 'sagemaker' \
    'unixodbc' && \
    mamba clean --all -f -y 

COPY add_r_kernel.sh /

RUN chmod +x /add_r_kernel.sh && \
    /add_r_kernel.sh && \
    chown ${NB_UID}:${NB_GID} /usr/local/lib/R -R

RUN usermod -a -G  sudo ${NB_USER} && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

WORKDIR $HOME
USER ${NB_UID}


CMD jupyter lab --ip 0.0.0.0 --port 8888 \
  --ServerApp.base_url="/jupyterlab/default" \
  --ServerApp.token='' \
  --ServerApp.allow_origin='*'
