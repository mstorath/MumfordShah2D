# Fixtures for the mumfordshah2d test suite.
# Ported by Claude Sonnet coding agent, Anthropic, 2026.

import numpy as np
import pytest


@pytest.fixture
def rng():
    return np.random.default_rng(seed=42)


@pytest.fixture
def step_signal_1d():
    """Clean step signal: 50 zeros followed by 50 ones."""
    return np.array([0.0] * 50 + [1.0] * 50)


@pytest.fixture
def piecewise_smooth_1d():
    """Piecewise-smooth signal with a single jump in the middle, sampled
    on [0, 1] with n=40 points: a smooth ramp on [0, 0.5), then a smooth
    parabola on [0.5, 1].
    """
    x = np.linspace(0.0, 1.0, 40)
    y = np.where(x < 0.5, 0.4 * x, 0.5 + 0.6 * (x - 0.5) ** 2)
    return y


@pytest.fixture
def step_image_2d():
    """64x64 grayscale image: left half = 0, right half = 1."""
    img = np.zeros((64, 64), dtype=np.float64)
    img[:, 32:] = 1.0
    return img


@pytest.fixture
def rgb_step_image():
    """64x64x3 RGB image: left half = (0,0,0), right half = (1,0.5,0.2)."""
    img = np.zeros((64, 64, 3), dtype=np.float64)
    img[:, 32:, 0] = 1.0
    img[:, 32:, 1] = 0.5
    img[:, 32:, 2] = 0.2
    return img
