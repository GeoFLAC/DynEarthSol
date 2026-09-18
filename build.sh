#!/bin/bash

# Set the number of dimensions
NDIMS=2 # or 3
# Set the GCC version
CXXVERSION=gcc-11 # clang-14 or gcc-8
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

# Pull the base image
docker pull $BASE_IMAGE
# Build the docker image
docker build --rm -t dynearthsol/$CXXVERSION . \
    --build-arg CXXVERSION=$CXXVERSION \
    --build-arg BASE_IMAGE=$BASE_IMAGE \
    --build-arg NDIMS=$NDIMS \
    --build-arg TZ=$HOST_TZ \

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Docker build succeeded."
    echo ""
    echo "To run the container execute:"
    echo "$ docker run -it --rm dynearthsol/$CXXVERSION"
else
    echo "Docker build failed. Cleaning up..."
    docker rmi dynearthsol/$CXXVERSION
fi

