"""Analysis-only spatial observables; none of these define a physical rate.

Array values are supplied by the caller. Lengths are lattice spacings. Connected
components use nearest-neighbor edges on the same periodic square as the model.
The spectral length is descriptive and need not be an intrinsic droplet radius.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np


def observe(material: Any, enzyme: Any, *, sticky: int = 2) -> dict[str, Any]:
    material = np.asarray(material)
    enzyme = np.asarray(enzyme)
    if material.ndim != 2 or material.shape[0] != material.shape[1]:
        raise ValueError("material must be a square two-dimensional field")
    if enzyme.shape != material.shape or not np.all((enzyme == 0) | (enzyme == 1)):
        raise ValueError("enzyme must be a matching binary field")
    side = material.shape[0]
    if side < 4:
        raise ValueError("side must be at least four")
    b = material == sticky
    total = side * side
    count = int(b.sum())
    visited = np.zeros_like(b)
    sizes: list[int] = []
    for row, col in np.argwhere(b):
        if visited[row, col]:
            continue
        visited[row, col] = True
        stack = [(int(row), int(col))]
        size = 0
        while stack:
            r, c = stack.pop()
            size += 1
            for rr, cc in (((r - 1) % side, c), ((r + 1) % side, c),
                           (r, (c - 1) % side), (r, (c + 1) % side)):
                if b[rr, cc] and not visited[rr, cc]:
                    visited[rr, cc] = True
                    stack.append((rr, cc))
        sizes.append(size)
    sizes.sort(reverse=True)
    oriented_bonds = sum(int(np.sum(b & np.roll(b, shift, axis)))
                         for axis in (0, 1) for shift in (-1, 1))
    neighbor_fraction = oriented_bonds / (4 * count) if count else None
    random_neighbor_fraction = (count - 1) / (total - 1) if count else None
    e_in = int(enzyme[b].sum())
    e_out = int(enzyme[~b].sum())
    pin = e_in / count if count else None
    pout = e_out / (total - count) if count < total else None
    enrichment = pin / pout if pin is not None and pout else None

    centered = b.astype(float) - count / total
    power = np.abs(np.fft.fft2(centered)) ** 2 / total
    k = 2 * math.pi * np.fft.fftfreq(side)
    kr = np.sqrt(k[:, None] ** 2 + k[None, :] ** 2)
    nonzero = kr > 0
    low = nonzero & (kr <= math.pi / 2 + 1e-12)
    spectral_mass = float(power[nonzero].sum())
    low_mass = float(power[low].sum())
    weighted_k = float((kr[low] * power[low]).sum())
    wavelength = 2 * math.pi * low_mass / weighted_k if weighted_k else None
    # The zero mode is excluded; all modes enter the Parseval check.
    parseval = float(power.sum())
    expected_parseval = count * (1 - count / total)
    if abs(parseval - expected_parseval) > 1e-9 * max(1, total):
        raise ArithmeticError("spectral normalization failed")
    return {
        "side": side, "n_B": count, "B_fraction": count / total,
        "n_E": int(enzyme.sum()), "cluster_sizes": sizes,
        "largest_cluster": sizes[0] if sizes else 0,
        "largest_B_fraction": sizes[0] / count if count else None,
        "mass_weighted_cluster_size": sum(s * s for s in sizes) / count if count else None,
        "clusters_ge_9": sum(s >= 9 for s in sizes),
        "B_mass_in_clusters_ge_9": sum(s for s in sizes if s >= 9),
        "neighbor_B_fraction": neighbor_fraction,
        "excess_neighbor_B_fraction": (neighbor_fraction - random_neighbor_fraction)
            if neighbor_fraction is not None else None,
        "E_on_B": e_in, "E_on_nonB": e_out,
        "E_density_B": pin, "E_density_nonB": pout, "E_enrichment": enrichment,
        "E_on_B_excess": e_in - int(enzyme.sum()) * count / total,
        "E_on_B_excess_fraction": (e_in / int(enzyme.sum()) - count / total)
            if int(enzyme.sum()) else 0.0,
        "spectral_wavelength_low_k": wavelength,
        "low_k_power_fraction": low_mass / spectral_mass if spectral_mass else None,
        "spectral_parseval_residual": parseval - expected_parseval,
    }


def fixtures() -> dict[str, Any]:
    """Deterministic sanity checks, separate from physical outcome evidence."""
    empty = np.zeros((8, 8), dtype=int)
    e = np.zeros_like(empty)
    assert observe(empty, e)["cluster_sizes"] == []
    full = np.full_like(empty, 2)
    assert observe(full, e)["cluster_sizes"] == [64]
    assert observe(full, e)["spectral_wavelength_low_k"] is None
    crossed = empty.copy()
    crossed[0, 0] = crossed[0, -1] = 2
    assert observe(crossed, e)["cluster_sizes"] == [2]
    isolated = empty.copy()
    isolated[0, 0] = isolated[4, 4] = 2
    assert observe(isolated, e)["cluster_sizes"] == [1, 1]
    block = empty.copy()
    block[:3, :3] = 2
    e[0, 0] = 1
    e[4, 4] = 1
    obs = observe(block, e)
    assert obs["cluster_sizes"] == [9]
    assert obs["clusters_ge_9"] == 1
    assert abs(obs["E_enrichment"] - 55 / 9) < 1e-12
    assert abs(obs["E_on_B_excess"] - (1 - 18 / 64)) < 1e-12
    translated = observe(np.roll(block, (3, 5), axis=(0, 1)),
                         np.roll(e, (3, 5), axis=(0, 1)))
    for key in obs:
        if isinstance(obs[key], float):
            assert math.isclose(obs[key], translated[key], abs_tol=1e-12)
        else:
            assert obs[key] == translated[key]
    return {"passed": True, "cases": 6, "translation_invariance": True}


if __name__ == "__main__":
    import json
    print(json.dumps(fixtures(), sort_keys=True))
