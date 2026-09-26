# This builds the image locally for the current platform

# Use the Ubuntu 22.04 release as the base image
ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}
# Set the maintainer
LABEL maintainer="chaseshyu@gmail.com"

# Set environment variable to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

ARG TZ=UTC
ARG CXXVERSION=gcc-11
# GOSPL=1 adds a conda `gospl` environment and gospl_extensions, and builds the
# 3D executable with use_gospl=1 (GoSPL coupling is 3D only).
ARG GOSPL=0
# Update package list and install necessary software
RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-program-options-dev \
    python3 \
    python3-pip \
    git \
    curl \
    ca-certificates \
    vim \
    sudo

RUN if [ "$CXXVERSION" = "gcc-8" ]; then \
      apt-get install -y software-properties-common \
      && add-apt-repository ppa:ubuntu-toolchain-r/test \
      && apt-get update \
      && apt-get install -y gcc-8 g++-8 \
      && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 \
      && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-8 60; \
    elif [ "$CXXVERSION" = "gcc-11" ]; then \
      apt-get install -y tzdata; \
    elif [ "$CXXVERSION" = "clang-14" ]; then \
      apt-get update && apt-get install -y \
      clang \
      libomp-dev \
      tzdata; \
    # elif [ "$CXXVERSION" = "clang-17" ]; then \
    #   apt-get install -y \
    #   wget tzdata \
    #   lsb-release gnupg \
    #   software-properties-common \
    #   cmake \
    #   && wget https://apt.llvm.org/llvm.sh && chmod +x llvm.sh && ./llvm.sh 17 \
    #   && add-apt-repository ppa:ubuntu-toolchain-r/test \
    #   && apt-get update && apt-get install -y \
    #   clang-17 \
    #   libomp-17-dev \
    #   && update-alternatives --install /usr/bin/clang clang /usr/bin/clang-17 100 \
    #   && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-17 100; \
    fi \
    && apt-get clean \
    # set time zone of system
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# Add and enable the default user
ARG USER=human
ARG USER_ID=1000
RUN adduser --disabled-password --gecos '' --uid $USER_ID $USER \
    && adduser $USER sudo \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers \
    && chown -R $USER:$USER /home/$USER

# Set $HOME
USER $USER
ENV HOME=/home/$USER
ENV USER=$USER

# pip install requirements
RUN pip install numpy scipy

# GoSPL stack: Miniforge, the gospl env (geodels channel, Python 3.11 to match
# the Makefile default), and gospl_extensions built against that env. Kept
# before the source copy so a code change does not rebuild this layer.
RUN if [ "$GOSPL" = "1" ]; then \
      curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o /tmp/miniforge.sh \
      && bash /tmp/miniforge.sh -b -p $HOME/miniforge3 && rm /tmp/miniforge.sh \
      && $HOME/miniforge3/bin/mamba create -y -n gospl -c geodels -c conda-forge gospl python=3.11 \
      && $HOME/miniforge3/bin/conda clean -afy \
      && git clone --depth 1 https://github.com/GeoFLAC/gospl_extensions.git $HOME/opt/gospl_extensions \
      && cd $HOME/opt/gospl_extensions/cpp_interface \
      && $HOME/miniforge3/bin/conda run -n gospl make install-local \
      && printf '. %s/miniforge3/etc/profile.d/conda.sh\nconda activate gospl\n' "$HOME" | sudo tee /etc/profile.d/gospl.sh > /dev/null \
      && printf '. %s/miniforge3/etc/profile.d/conda.sh\nconda activate gospl\n' "$HOME" >> ~/.bashrc; \
    fi

# Clone the repository
COPY --chown=$USER:$USER . $HOME/DynEarthSol

# make dynearthsol2d
WORKDIR $HOME/DynEarthSol
ARG NDIMS=2
RUN make ndims=$NDIMS cleanall \
    && if [ "$GOSPL" = "1" ]; then \
      make ndims=3 use_gospl=1 -j4 \
        CONDA_ENV_PATH=$HOME/miniforge3/envs/gospl PYTHON_VERSION=3.11; \
    elif [ "$CXXVERSION" = "clang-14" ]; then \
      make ndims=$NDIMS -j4 CXX=clang++;\
    else \
      make ndims=$NDIMS -j4; \
    fi
# Default use 8 threads
ENV OMP_NUM_THREADS=8

# Set alias
RUN echo '#!/bin/bash\nexport OMP_NUM_THREADS=$1 && echo OMP_NUM_THREADS = $OMP_NUM_THREADS' \
    | sudo tee /usr/local/bin/set_omp_threads > /dev/null \
    && sudo chmod +x /usr/local/bin/set_omp_threads \
    && echo 'function omp() { . /usr/local/bin/set_omp_threads $1; }' >> ~/.bashrc

CMD ["bash"]