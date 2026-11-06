#!/bin/bash
# Install what a DES build needs beyond the compiler in a Debian-based gcc:N image
# (run as root): cmake, which configures the mmg submodule.
set -euo pipefail

# gcc:8-11 are Debian buster/bullseye: buster's packages are on archive.debian.org,
# and bullseye-security serves no pool files. A CI build needs neither suite.
if [ -f /etc/apt/sources.list ]; then
    sed -i -e '/security/d' -e '/-updates/d' /etc/apt/sources.list
    sed -i '/buster/s|deb.debian.org|archive.debian.org|' /etc/apt/sources.list
fi
apt-get -o Acquire::Check-Valid-Until=false update
apt-get install -y --no-install-recommends cmake
