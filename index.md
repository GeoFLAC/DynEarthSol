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
  <a class="video-thumb" href="https://youtu.be/0AuTI7MZ05M" data-video-id="0AuTI7MZ05M" target="_blank" rel="noopener">
    <img src="https://img.youtube.com/vi/0AuTI7MZ05M/hqdefault.jpg" alt="Video 0AuTI7MZ05M thumbnail">
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

<p style="text-align: center;">
  <a href="https://www.youtube.com/playlist?list=PLLxHSgtMeGthxjsfFTv4aPTzCfcOqNCYD">More simulation videos on YouTube &rarr;</a>
</p>

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

# Start here

{% include feature_row feature_row=site.data.features.start_row %}

# Quick start

Build the 2D executable and run a worked example. The
[user manual](https://geoflac.github.io/des3d/docs/usage) carries the complete
build guide: every build option, the optional libraries, and the platform notes.

**1. Install the dependencies.** Pick the line for your platform; the build
finds everything on its own, with no paths to edit.

```sh
brew install boost libomp                               # macOS (Homebrew)
sudo apt install g++ make libboost-program-options-dev  # Debian / Ubuntu
sudo dnf install gcc-c++ make boost-devel               # Fedora / RHEL
```

**2. Get the source**, including the bundled submodules.

```sh
git clone --recurse-submodules https://github.com/GeoFLAC/DynEarthSol.git
cd DynEarthSol
```

**3. Build.** `make` on its own builds the 3D executable, `dynearthsol3d`.
Pass `ndims=2` for the 2D executable used below.

```sh
make ndims=2
```

**4. Run an example.**

```sh
cd examples
../dynearthsol2d ./core-complex.cfg
```

More example configurations live in `examples/`, and `examples/defaults.cfg`
documents every input parameter. `dynearthsol2d -h` lists them too. To write
your own input file without editing one by hand, use the web-based
[input generator](https://geoflac.github.io/des-inputgen/), which builds a
ready-to-run `.cfg` file from a form. If a build
picks up an unexpected library, `make config` prints the compiler and every
dependency path it resolved to, which is also the most useful thing to include
in a bug report.

# Citing DES3D

Please cite the DES3D v2.0 paper (preprint):

> Shyu, C. J., Lee, S., Ding, X., Keum, J., Tan, E., Lavier, L. L., & Choi, E.
> (2026). DynEarthSol v2.0: An efficient explicit Lagrangian solver for
> geodynamics, surface processes, and earthquake-cycle dynamics. *EGUsphere*
> [preprint]. [doi:10.5194/egusphere-2026-2922](https://doi.org/10.5194/egusphere-2026-2922)

and the original method paper:

> Choi, E., Tan, E., Lavier, L. L., & Calo, V. M. (2013). DynEarthSol2D: An
> efficient unstructured finite element method to study long-term tectonic
> deformation. *Journal of Geophysical Research: Solid Earth*, 118(5),
> 2429–2444. [doi:10.1002/jgrb.50148](https://doi.org/10.1002/jgrb.50148)

To cite the software itself, use
[doi:10.5281/zenodo.20293557](https://doi.org/10.5281/zenodo.20293557). This
DOI covers every archived release, both the software (`v2.x.x`) and the
benchmark datasets (`benchmarks-v2.x.x`), and it opens the list of all
versions. To cite the exact version you used, pick the matching software
release from that list, or take its DOI from the
[releases page](https://github.com/GeoFLAC/DynEarthSol/releases) or from the
`CITATION.cff` file shipped with that release.


