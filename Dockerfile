# Builds the DynEarthSol image for the current platform.
#
# Stages, most shared first: base (compilers), env-0 (plain) or env-1 (GoSPL),
# then the user, timezone and source. Changing the uid, timezone or source
# rebuilds only the last stage, not the conda environment.

# GOSPL=1: add the gospl conda env and build the 3D GoSPL executable.
ARG GOSPL=0

FROM ubuntu:22.04 AS base
ENV DEBIAN_FRONTEND=noninteractive
ARG CXXVERSION=gcc-11
# gcc-11 is the Ubuntu 22.04 default; clang-14 adds clang and its OpenMP.
# apt downloads go to cache mounts, kept between builds but not in the image.
# docker-clean is removed because it would empty that cache.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean \
    && apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-program-options-dev \
    python3 \
    python3-pip \
    git \
    curl \
    ca-certificates \
    tzdata \
    vim \
    sudo \
    && if [ "$CXXVERSION" = "clang-14" ]; then \
      apt-get install -y clang libomp-dev; \
    fi

# Plain image: HDF5 for hdf5=1, and numpy/scipy for the Python utilities.
FROM base AS env-0
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    --mount=type=cache,target=/root/.cache/pip \
    apt-get update && apt-get install -y libhdf5-dev pkg-config \
    && pip install numpy scipy

# GoSPL image: Miniforge, the gospl env (Python 3.11, the Makefile default)
# and gospl_extensions, in /opt, writable by the `conda` group.
# - HDF5 comes from the env, not apt: GoSPL already loads the env's libhdf5
#   into the DES process.
# - Versions are pinned to the last set tested together (2026-11-04);
#   override them with --build-arg.
# - No cache mount for conda: it then copies instead of hard-linking, and
#   the image grew by 110 MB.
# - curl --retry-all-errors also retries a dropped download. 5 is a guess.
# - CONDA_OVERRIDE_ARCHSPEC=x86_64: build for any amd64 CPU, not the build
#   machine's. An AVX-512 env from the layer cache crashed (SIGILL) elsewhere.
FROM base AS env-1
ARG MINIFORGE_VERSION=26.9.4-1
ARG GOSPL_VERSION=2026.7.14
ARG GOSPL_EXT_COMMIT=ea0c01de3e9eb80bbf224ea644d21c72d073a390
RUN curl -fsSL --retry 5 --retry-all-errors \
      "https://github.com/conda-forge/miniforge/releases/download/$MINIFORGE_VERSION/Miniforge3-$MINIFORGE_VERSION-Linux-$(uname -m).sh" -o /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/miniforge3 && rm /tmp/miniforge.sh \
    && if [ "$(uname -m)" = x86_64 ]; then export CONDA_OVERRIDE_ARCHSPEC=x86_64; fi \
    && /opt/miniforge3/bin/mamba create -y -n gospl -c geodels -c conda-forge gospl=$GOSPL_VERSION python=3.11 \
    && /opt/miniforge3/bin/conda clean -afy \
    && git init -q /opt/gospl_extensions \
    && cd /opt/gospl_extensions \
    && git remote add origin https://github.com/GeoFLAC/gospl_extensions.git \
    && git fetch -q --depth 1 origin $GOSPL_EXT_COMMIT \
    && git checkout -q FETCH_HEAD \
    && cd cpp_interface \
    && /opt/miniforge3/bin/conda run -n gospl make install-local \
    && groupadd conda \
    && chgrp -R conda /opt/miniforge3 /opt/gospl_extensions \
    && chmod -R g+rwX /opt/miniforge3 /opt/gospl_extensions \
    && printf '. /opt/miniforge3/etc/profile.d/conda.sh\nconda activate gospl\nexport HDF5_INCLUDE_DIR=$CONDA_PREFIX/include HDF5_LIB_DIR=$CONDA_PREFIX/lib\n' \
      | tee /etc/profile.d/gospl.sh >> /etc/bash.bashrc

FROM env-${GOSPL}
LABEL maintainer="chaseshyu@gmail.com"
ARG GOSPL

# The user. In the GoSPL image it joins `conda`, and ~/miniforge3 and
# ~/opt/gospl_extensions link to /opt, where the Makefile looks.
# useradd -l: skip lastlog, a file sized by the uid (uid 5000000 added 1.6 GB).
ARG USER=human
ARG USER_ID=1000
RUN useradd -l -m -s /bin/bash -u $USER_ID -G sudo $USER \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers \
    && if [ "$GOSPL" = "1" ]; then \
      usermod -aG conda $USER \
      && mkdir -p /home/$USER/opt \
      && ln -s /opt/miniforge3 /home/$USER/miniforge3 \
      && ln -s /opt/gospl_extensions /home/$USER/opt/gospl_extensions; \
    fi \
    && chown -R $USER:$USER /home/$USER

# Timezone
ARG TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

USER $USER
ENV HOME=/home/$USER
ENV USER=$USER

# The source; .dockerignore leaves out the host's build output.
COPY --chown=$USER:$USER . $HOME/DynEarthSol

WORKDIR $HOME/DynEarthSol
# Fetch the submodules now, so a build in the container needs no network.
# knn-bvh is left out: only openacc=1 uses it, and the image has no NVHPC.
RUN git submodule update --init --depth 1 mmg nanoflann
ARG CXXVERSION=gcc-11
ARG NDIMS=2
RUN if [ "$CXXVERSION" = "clang-14" ]; then cxx=CXX=clang++; fi \
    && make ndims=$NDIMS cleanall \
    && if [ "$GOSPL" = "1" ]; then \
      make ndims=3 usegospl=1 usemmg=1 -j$(nproc) $cxx \
        CONDA_ENV_PATH=$HOME/miniforge3/envs/gospl PYTHON_VERSION=3.11; \
    else \
      make ndims=$NDIMS -j$(nproc) $cxx; \
    fi
# OMP_NUM_THREADS is unset, so OpenMP uses every CPU; `omp N` sets it.
RUN echo '#!/bin/bash\nexport OMP_NUM_THREADS=$1 && echo OMP_NUM_THREADS = $OMP_NUM_THREADS' \
    | sudo tee /usr/local/bin/set_omp_threads > /dev/null \
    && sudo chmod +x /usr/local/bin/set_omp_threads \
    && echo 'function omp() { . /usr/local/bin/set_omp_threads $1; }' >> ~/.bashrc

CMD ["bash"]
