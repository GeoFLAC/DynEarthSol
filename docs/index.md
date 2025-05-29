---
title: "DynEarthSol"
layout: splash
permalink: /
header:
  overlay_image: /assets/images/banner_lowres.png
tagline: "Dynamic Earth Solver, a standard finite element transplant of geoflac for unstructured meshes with P1 elements in 2D and 3D."
---

---
layout: archive
title: "Main Features"
<!-- permalink: /post-archive-feature-rows/ -->
<!-- author_profile: true -->
feature_row:
  - image_path: /assets/images/Globe.svg
    alt: "Globe"
    title: "Explicit dynamic relaxation"
    description: "DynEarthSol finds static or quasi-static equilibrium through dynamic relaxation. Combined with explicit time marching, this feature significantly improves stability and allows for large deformations."
  - image_path: /assets/images/Globe.svg
    alt: "Globe"
    title: "Updated Lagrangian scheme"
    description: "Can simulate large deformations within the framework of small-strain kinematics and keeps track of free boundary deformation."
  - image_path: /assets/images/Globe.svg
    alt: "Globe"
    title: "Various material models"
    description: "Provides linear elasticity, poro-elasticity, Mohr-Coulomb plasticity, and Maxwell viscoelasticity."
<!-- {% for post in site.posts limit: 5 %}
  {% include archive-single.html %}
{% endfor %} -->
{% include feature_row id="intro" type="center" %}
{% include feature_row %}
<!-- {% include feature_row id="main_features" %} -->
---

# Get Started

For full instruction on how to install DynEarthSol, please visit the [online user manual](https://geoflac.github.io/des3d/docs/usage).

## Step 1: Check prerequisites

- A recent C++ compiler supporting C++11 (g++ 4.4+)
- Boost::Program_options library (1.42+)
- Python 2.6+ or 3.2+ with Numpy

## Step 2: Get DynEarthSol

```sh
git clone https://github.com/GeoFLAC/DynEarthSol
```

## Step 3: Build

```sh
cd DynEarthSol
make
```

## Step 4: Run

```sh
cd examples
../dynearthsol2d ./core-complex.cfg
```


