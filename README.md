# MumfordShah2D — Edge-preserving image restoration via the Mumford-Shah model

[![PyPI](https://img.shields.io/pypi/v/mumfordshah2d.svg)](https://pypi.org/project/mumfordshah2d/)
[![Python](https://img.shields.io/pypi/pyversions/mumfordshah2d.svg)](https://pypi.org/project/mumfordshah2d/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![CI](https://github.com/mstorath/MumfordShah2D/actions/workflows/ci.yml/badge.svg)](https://github.com/mstorath/MumfordShah2D/actions/workflows/ci.yml)
[![MATLAB](https://img.shields.io/badge/MATLAB-supported-orange.svg)](#matlab)

Edge-preserving image restoration via the Mumford-Shah model — denoising, deconvolution, and inpainting of vector-valued images with linear complexity in the number of channels and no a-priori label discretisation.

![Mumford-Shah denoising example](/Docs/img_salt_pepper.png)

## Paper

> K. Hohm, M. Storath, A. Weinmann.
> [*An algorithmic framework for Mumford-Shah regularization of inverse problems in imaging.*](http://bigwww.epfl.ch/publications/hohm1501.pdf)
> Inverse Problems 31(11), 115011, 2015.

## Quickstart

### Python (Rust core)

```bash
pip install mumfordshah2d
```

```python
import numpy as np
from mumfordshah2d import min_l2_mum_2d

# noisy grayscale image
f = np.random.randn(64, 64) + your_image
u = min_l2_mum_2d(f, gamma=0.5, alpha=1.0)
```

The Python package wraps a Rust extension built with [PyO3](https://pyo3.rs)
and [maturin](https://maturin.rs); algorithm crate lives under `src/`,
demos under [`demos_python/`](demos_python/).
See [`README_PYTHON.md`](README_PYTHON.md) for the full Python API,
including soft / hard thresholding utilities and prox handles.

### MATLAB

The original MATLAB / Java reference implementation is in this same repository:

1. Run `setPath.m` to add the necessary folders to the MATLAB path.
2. For best performance, increase Java heap space in the MATLAB preferences (MATLAB → General → Java Heap Memory).
3. Run a demo from the `Demos/` folder.

## Application examples

### Edge-preserving smoothing of vector-valued images

- Supports smoothing of vector-valued images (e.g. multispectral, feature images)
- Linear complexity in the number of channels
- No discretisation of colour space required

(See hero image above for the salt-and-pepper denoising case.)

### Regularization for deconvolution

![Deconvolution](/Docs/img_deconv.png)

### Inpainting

![Inpainting](/Docs/img_inpainting.png)

## How to cite

If you use this software, please cite the paper above. GitHub's "Cite this repository" button on the repo page reads the `version` and `date-released` fields from [`CITATION.cff`](CITATION.cff) and renders BibTeX/APA.

## See also

Sibling projects from the same research program on variational methods for signal and image processing:

- [Pottslab](https://github.com/mstorath/Pottslab) — multilabel image segmentation via the Potts / piecewise-constant Mumford-Shah model
- [L1TV](https://github.com/mstorath/L1TV) — exact L1-TV regularisation of real- or circle-valued signals
- [CSSD](https://github.com/mstorath/CSSD) — cubic smoothing splines for signals with discontinuities
- [CircleMedianFilter](https://github.com/mstorath/CircleMedianFilter) — fast median filtering for phase or orientation data
- [DCEBE](https://github.com/mstorath/DCEBE) — bolus arrival time estimation for DCE-MRI signals

## License

Released under the MIT License. See [LICENSE](LICENSE).

---

### Project history

The Python / Rust port of this codebase was generated from the original MATLAB / Java reference by a Claude coding agent in 2026.
See [`PORTED_BY.md`](PORTED_BY.md) for full attribution and the porting plan.
