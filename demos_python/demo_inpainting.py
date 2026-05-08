"""Mumford-Shah inpainting demo — Python port of Demos/demoInpainting.m.

Reproduces the MATLAB demo behaviour: load a colour image, corrupt it
with Gaussian noise and a random binary mask (60% missing pixels), then
restore via the 2-D L2 Mumford-Shah model with a weighted-L2 data
fidelity prox where missing pixels carry weight 0.

Run from the repo root::

    python demos_python/demo_inpainting.py [--smoke]

``--smoke`` reduces image size and iterations for fast CI smoke tests.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None  # demo runs headless without matplotlib; just skips plotting

from mumfordshah2d import min_l2_mum_2d
from mumfordshah2d.prox import make_prox_l2w


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IMAGE = REPO_ROOT / "Demos" / "fruitsColor.jpg"


def load_image(path: Path, smoke: bool) -> np.ndarray:
    """Load an image as float64 ``(rows, cols, channels)`` in [0, 1].

    Falls back to a synthetic patchy gradient if PIL/imageio aren't
    available or the file is missing — keeps the demo runnable in
    minimal environments and CI smoke contexts.
    """
    try:
        from PIL import Image  # local import: matplotlib draws PIL too

        img = np.asarray(Image.open(path).convert("RGB"), dtype=np.float64) / 255.0
    except Exception:
        # Synthetic fallback: 64×64 colour patches with sharp edges, useful
        # for both smoke testing and dependency-free runs.
        h = w = 32 if smoke else 96
        img = np.zeros((h, w, 3), dtype=np.float64)
        img[: h // 2, : w // 2] = [0.85, 0.20, 0.20]
        img[: h // 2, w // 2 :] = [0.20, 0.85, 0.20]
        img[h // 2 :, : w // 2] = [0.20, 0.20, 0.85]
        img[h // 2 :, w // 2 :] = [0.85, 0.85, 0.20]
        return img

    if smoke:
        # Subsample for fast CI runs.
        img = img[::4, ::4]
    return img


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="fast smoke run")
    parser.add_argument("--missing", type=float, default=0.6, help="fraction of pixels masked out")
    parser.add_argument("--sigma", type=float, default=0.02, help="Gaussian noise std")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--image", type=Path, default=DEFAULT_IMAGE)
    parser.add_argument("--save", type=Path, default=None, help="optional output path for the figure")
    args = parser.parse_args(argv)

    img = load_image(args.image, args.smoke)
    rng = np.random.default_rng(args.seed)
    noisy = np.clip(img + args.sigma * rng.standard_normal(img.shape), 0.0, 1.0)

    # MATLAB: weights = rand(m,n) > missingFraction; imgNoisy(~mask3) = 0
    weights = (rng.random(img.shape[:2]) > args.missing).astype(np.float64)
    corrupted = noisy * weights[:, :, np.newaxis]

    # MATLAB demo parameters: s = 0.07, alpha = 0.8, gamma = s^2 * alpha
    s = 0.07
    alpha = 0.8
    gamma = s * s * alpha
    max_iter = 200 if args.smoke else 5000
    tol = 0.05 if args.smoke else 0.01

    # Construct the prox over the *corrupted* image: weights=0 at missing
    # pixels means "ignore the data, take the consensus", i.e. the inpainted
    # value comes entirely from neighbouring known regions via the MS prior.
    prox = make_prox_l2w(corrupted, weights)

    print(f"Image shape: {img.shape}, missing pixels: {(weights == 0).mean():.1%}")
    print(f"Solving L2-Mumford-Shah: gamma={gamma:.4g}, alpha={alpha}, max_iter={max_iter}")

    restored = min_l2_mum_2d(
        corrupted,
        gamma=gamma,
        alpha=alpha,
        tol=tol,
        max_iter=max_iter,
        isotropic=2,  # knight-move (matches MATLAB demo)
        prox=prox,
        verbose=not args.smoke,
    )

    psnr_corrupted = -10 * np.log10(np.mean((noisy - corrupted) ** 2) + 1e-12)
    psnr_restored = -10 * np.log10(np.mean((img - restored) ** 2) + 1e-12)
    print(f"PSNR corrupted vs noiseless: {psnr_corrupted:.2f} dB")
    print(f"PSNR restored  vs original : {psnr_restored:.2f} dB")

    if plt is None or args.smoke:
        return 0

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, im, title in zip(
        axes,
        [img, corrupted, np.clip(restored, 0, 1)],
        ["Original", "Noisy + 60% missing", "Mumford-Shah restoration"],
    ):
        ax.imshow(im)
        ax.set_title(title)
        ax.axis("off")
    fig.tight_layout()
    if args.save is not None:
        fig.savefig(args.save, dpi=120)
        print(f"Saved figure to {args.save}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
