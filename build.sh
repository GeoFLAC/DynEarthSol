#!/bin/bash

# Each setting can also be passed in the environment.
NDIMS=${NDIMS:-2}                # or 3
CXXVERSION=${CXXVERSION:-gcc-11} # or clang-14
GOSPL=${GOSPL:-0}                # 1: GoSPL coupling, 3D only
HOST_TZ=$(readlink /etc/localtime | sed 's|.*/zoneinfo/||')
# Your uid, so the container can write ./runs. Under sudo, 1000: uid 0 exists.
USER_ID=${USER_ID:-$(id -u)}
if [ "$USER_ID" = "0" ]; then
    USER_ID=1000
fi

# The Dockerfile's cache mounts need BuildKit.
if ! docker buildx version > /dev/null 2>&1; then
    echo "build.sh needs Docker with BuildKit (the docker buildx plugin)."
    exit 1
fi
# The image fetches the submodules through .git, so a worktree cannot work.
if [ ! -d .git ]; then
    echo "build.sh needs a regular clone, run from its top folder: .git is not a directory here."
    exit 1
fi

if [ "$CXXVERSION" != "gcc-11" ] && [ "$CXXVERSION" != "clang-14" ]; then
    echo "CXXVERSION must be gcc-11 or clang-14, not $CXXVERSION."
    exit 1
fi

VARIANT=$CXXVERSION
if [ "$GOSPL" = "1" ]; then
    NDIMS=3
    VARIANT=$VARIANT-gospl
fi
TAG=dynearthsol/$VARIANT

# Layers CI built on master (.github/workflows/docker.yml). On a miss, the
# build warns and builds them here.
case $(uname -m) in
    x86_64) ARCH=amd64 ;;
    arm64|aarch64) ARCH=arm64 ;;
    *) ARCH=$(uname -m) ;;
esac
CACHE=ghcr.io/geoflac/dynearthsol/docker-cache:$VARIANT-$ARCH

docker pull ubuntu:22.04
docker build --rm -t $TAG . \
    --cache-from type=registry,ref=$CACHE \
    --build-arg CXXVERSION=$CXXVERSION \
    --build-arg NDIMS=$NDIMS \
    --build-arg GOSPL=$GOSPL \
    --build-arg TZ=$HOST_TZ \
    --build-arg USER_ID=$USER_ID

# A failed build keeps the previous image under the tag.
if [ $? -ne 0 ]; then
    echo "Docker build failed; $TAG is unchanged."
    exit 1
fi
echo "Docker build succeeded."
echo ""
echo "To run the container execute:"
echo "$ mkdir -p runs && docker run -it --rm --network none -v \"\$PWD/runs\":/home/human/DynEarthSol/runs $TAG"
echo "Only results in ~/DynEarthSol/runs (./runs here) are kept after exit."

