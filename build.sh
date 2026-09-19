#!/bin/bash

# Set the number of dimensions
NDIMS=${NDIMS:-2} # or 3; every knob here can also be set in the environment
# Set the GCC version
CXXVERSION=${CXXVERSION:-gcc-11} # clang-14 or gcc-8
# Set to 1 for GoSPL surface-process coupling: adds a conda gospl environment
# and gospl_extensions, forces NDIMS=3, and tags the image <CXXVERSION>-gospl
GOSPL=${GOSPL:-0}
# Set the timezone of the host machine
HOST_TZ=$(readlink /etc/localtime | sed 's|.*/zoneinfo/||')

# based on gcc version to choose the image
if [ "$CXXVERSION" = "gcc-8" ]; then
    BASE_IMAGE=ubuntu:20.04
elif [ "$CXXVERSION" = "gcc-11" ]; then
    BASE_IMAGE=ubuntu:22.04
elif [ "$CXXVERSION" = "clang-14" ]; then
    BASE_IMAGE=ubuntu:22.04
else
    echo "CXX version not supported"
    exit 1
fi

TAG=dynearthsol/$CXXVERSION
if [ "$GOSPL" = "1" ]; then
    NDIMS=3
    TAG=$TAG-gospl
fi

# Pull the base image
docker pull $BASE_IMAGE
# Build the docker image
docker build --rm -t $TAG . \
    --build-arg CXXVERSION=$CXXVERSION \
    --build-arg BASE_IMAGE=$BASE_IMAGE \
    --build-arg NDIMS=$NDIMS \
    --build-arg GOSPL=$GOSPL \
    --build-arg TZ=$HOST_TZ \

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Docker build succeeded."
    echo ""
    echo "To run the container execute:"
    echo "$ docker run -it --rm $TAG"
else
    echo "Docker build failed. Cleaning up..."
    docker rmi $TAG
fi

