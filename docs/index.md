---
title: "DynEarthSol"
layout: splash
permalink: /
header:
  overlay_image: /assets/images/banner_lowres.png
  tagline: "Dynamic Earth Solver, a standard finite element transplant of geoflac for unstructured meshes with P1 elements in 2D and 3D."
---

{% include feature_row feature_row=site.data.features.feature_row %}
<img src="assets/images/nsf-logo.png" alt="NSF Logo" width="32"/>Supported by the National Science Foundation Award 2104002.

<div class="video-grid">
  <a class="video-thumb" href="https://youtu.be/jwDE1u7aFjY" data-video-id="jwDE1u7aFjY" target="_blank" rel="noopener">
    <img src="https://img.youtube.com/vi/jwDE1u7aFjY/hqdefault.jpg" alt="Video jwDE1u7aFjY thumbnail">
  </a>
  <a class="video-thumb" href="https://youtu.be/vfeiY5WuV9c" data-video-id="vfeiY5WuV9c" target="_blank" rel="noopener">
    <img src="https://img.youtube.com/vi/vfeiY5WuV9c/hqdefault.jpg" alt="Video vfeiY5WuV9c thumbnail">
  </a>
  <a class="video-thumb" href="https://youtu.be/zr-4HIg7_14" data-video-id="zr-4HIg7_14" target="_blank" rel="noopener">
    <img src="https://img.youtube.com/vi/zr-4HIg7_14/hqdefault.jpg" alt="Video zr-4HIg7_14 thumbnail">
  </a>
  <a class="video-thumb" href="https://youtu.be/_VZUCwzfxPk" data-video-id="_VZUCwzfxPk" target="_blank" rel="noopener">
    <img src="https://img.youtube.com/vi/_VZUCwzfxPk/hqdefault.jpg" alt="Video _VZUCwzfxPk thumbnail">
  </a>
</div>

<script>
  document.addEventListener('DOMContentLoaded', function () {
    var thumbs = document.querySelectorAll('.video-grid .video-thumb');
    thumbs.forEach(function (link) {
      link.addEventListener('click', function (e) {
        e.preventDefault();
        var vid = link.getAttribute('data-video-id');
        if (!vid) return;
        var embedUrl = 'https://www.youtube-nocookie.com/embed/' + vid + '?autoplay=1&rel=0';
        var container = document.createElement('div');
        container.className = 'responsive-video-container';
        var iframe = document.createElement('iframe');
        iframe.src = embedUrl;
        iframe.setAttribute('frameborder', '0');
        iframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share');
        iframe.setAttribute('allowfullscreen', 'true');
        container.appendChild(iframe);
        link.replaceWith(container);
      }, { once: true });
    });
  });
  </script>

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


