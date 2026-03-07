#!/usr/bin/env bash
# install-openmp.sh — Download, build, and install LLVM OpenMP for macOS.
# Usage: bash install-openmp.sh <LLVM_VERSION>
# Example: bash install-openmp.sh 19.1.7

set -e

# --- Argument check ---
if [ $# -ne 1 ]; then
    echo "Usage: $(basename "$0") <LLVM_VERSION>" >&2
    echo "Example: $(basename "$0") 19.1.7" >&2
    exit 1
fi

VERSION="$1"

# --- cmake check ---
if ! command -v cmake >/dev/null 2>&1; then
    echo "Error: cmake not found. Please install cmake first." >&2
    exit 1
fi

# --- Absolute path to this script's directory ---
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# --- Architecture detection ---
ARCH=$(uname -m)
case "$ARCH" in
    arm64)   ARCH=arm64 ;;
    x86_64)  ARCH=x86_64 ;;
    *)
        echo "Warning: Unknown architecture '$ARCH', falling back to x86_64." >&2
        ARCH=x86_64
        ;;
esac

# --- Paths ---
INSTALL_PREFIX="$SCRIPT_DIR/openmp-install"
CMAKE_DIR="$SCRIPT_DIR/cmake-${VERSION}.src"
OPENMP_SRC_DIR="$SCRIPT_DIR/openmp-${VERSION}.src"

BASE_URL="https://github.com/llvm/llvm-project/releases/download/llvmorg-${VERSION}"
CMAKE_TARBALL="$SCRIPT_DIR/cmake-${VERSION}.src.tar.xz"
OPENMP_TARBALL="$SCRIPT_DIR/openmp-${VERSION}.src.tar.xz"

# --- Download cmake tarball ---
if [ ! -f "$CMAKE_TARBALL" ]; then
    echo "Downloading cmake-${VERSION}.src.tar.xz ..."
    curl -L "${BASE_URL}/cmake-${VERSION}.src.tar.xz" -o "$CMAKE_TARBALL"
else
    echo "cmake-${VERSION}.src.tar.xz already downloaded, skipping."
fi

# --- Download openmp tarball ---
if [ ! -f "$OPENMP_TARBALL" ]; then
    echo "Downloading openmp-${VERSION}.src.tar.xz ..."
    curl -L "${BASE_URL}/openmp-${VERSION}.src.tar.xz" -o "$OPENMP_TARBALL"
else
    echo "openmp-${VERSION}.src.tar.xz already downloaded, skipping."
fi

# --- Extract cmake tarball ---
if [ ! -d "$CMAKE_DIR" ]; then
    echo "Extracting cmake-${VERSION}.src.tar.xz ..."
    tar xf "$CMAKE_TARBALL" -C "$SCRIPT_DIR"
else
    echo "cmake-${VERSION}.src already extracted, skipping."
fi

# --- Extract openmp tarball ---
if [ ! -d "$OPENMP_SRC_DIR" ]; then
    echo "Extracting openmp-${VERSION}.src.tar.xz ..."
    tar xf "$OPENMP_TARBALL" -C "$SCRIPT_DIR"
else
    echo "openmp-${VERSION}.src already extracted, skipping."
fi

# --- Build ---
mkdir -p "${OPENMP_SRC_DIR}/build"
cd "${OPENMP_SRC_DIR}/build"

# Clear compiler/linker env vars that can leak in from conda or Make and
# cause CMake's compiler test to fail (e.g. -lomp injected before libomp exists).
unset LDFLAGS CFLAGS CXXFLAGS CPPFLAGS

cmake \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_MODULE_PATH="${CMAKE_DIR}/Modules" \
    -DCMAKE_OSX_ARCHITECTURES="$ARCH" \
    -DLIBOMP_INSTALL_ALIASES=OFF \
    ..

make -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
make install

echo "OpenMP ${VERSION} installed to external/openmp-install"
