---
title: "Quickstart"
layout: single
permalink: /quickstart/
---

For full instruction on how to install DynEarthSol, please visit the [online user manual](https://geoflac.github.io/des3d/docs/usage).

### Step 1: Check prerequisites

- A recent C++ compiler supporting C++11 (g++ 4.4+)
- Boost::Program_options library (1.42+)
- Python 2.6+ or 3.2+ with Numpy

### Step 2: Get DynEarthSol

```sh
git clone https://github.com/GeoFLAC/DynEarthSol
```

### Step 3: Build

```sh
cd DynEarthSol
make
```

### Step 4: Run

```sh
cd examples
../dynearthsol2d ./core-complex.cfg
```

