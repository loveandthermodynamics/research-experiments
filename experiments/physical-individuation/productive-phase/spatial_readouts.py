"""Observation-only morphology and template association; never natural rates."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("frozen_morphology", ROOT / "dependencies/morphology.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def dense_core(material, minimum=9):
    array = np.asarray(material)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("square grid required")
    side = array.shape[0]
    b, visited, core = array == 2, np.zeros_like(array, dtype=bool), np.zeros_like(array, dtype=bool)
    for rr, cc in np.argwhere(b):
        if visited[rr, cc]:
            continue
        visited[rr, cc] = True
        stack, component = [(int(rr), int(cc))], []
        while stack:
            r, c = stack.pop()
            component.append((r, c))
            for i, j in (((r-1) % side, c), ((r+1) % side, c),
                         (r, (c-1) % side), (r, (c+1) % side)):
                if b[i, j] and not visited[i, j]:
                    visited[i, j] = True
                    stack.append((i, j))
        if len(component) >= minimum:
            for i, j in component:
                core[i, j] = True
    return core


def observe(material, templates):
    material, templates = np.asarray(material), np.asarray(templates)
    if material.shape != templates.shape or not np.isin(templates, (0, 1, 2)).all():
        raise ValueError("matching physical template layer required")
    # The old binary enzyme argument is only a Y-occupancy readout here.
    result = module.observe(material, templates == 2)
    for key in list(result):
        if key.startswith("E_"):
            result["Y_"+key[2:]] = result.pop(key)
    result["n_Y"] = result.pop("n_E")
    result["n_X"] = int((templates == 1).sum())
    result["n_templates"] = result["n_X"] + result["n_Y"]
    result["template_density"] = result["n_templates"] / material.size
    core = dense_core(material)
    assert int(core.sum()) == result["B_mass_in_clusters_ge_9"]
    near = core.copy()
    for axis in (0, 1):
        for shift in (-1, 1):
            near |= np.roll(core, shift, axis)
    present = templates != 0
    result["templates_near_dense_B"] = int((present & near).sum())
    result["template_density_near_dense_B"] = result["templates_near_dense_B"] / material.size
    result["X_near_dense_B"] = int(((templates == 1) & near).sum())
    result["Y_near_dense_B"] = int(((templates == 2) & near).sum())
    result["dense_neighborhood_sites"] = int(near.sum())
    result["large_B_mass_fraction"] = (result["B_mass_in_clusters_ge_9"] / result["n_B"]
                                        if result["n_B"] else 0.)
    return result


def fixtures():
    a = np.zeros((8, 8), dtype=int)
    t = np.zeros_like(a)
    a[:3, :3] = 2
    t[1, 1] = 1
    t[3, 1] = 2
    t[6, 6] = 1
    result = observe(a, t)
    assert result["templates_near_dense_B"] == 2
    assert result["X_near_dense_B"] == result["Y_near_dense_B"] == 1
    assert result["n_templates"] == 3 and result["large_B_mass_fraction"] == 1
    moved = observe(np.roll(a, (4, 5), (0, 1)), np.roll(t, (4, 5), (0, 1)))
    for key in result:
        if isinstance(result[key], float):
            assert abs(result[key]-moved[key]) < 1e-12
        else:
            assert result[key] == moved[key]
    none = observe(np.zeros_like(a), t)
    assert none["templates_near_dense_B"] == 0
    assert none["large_B_mass_fraction"] == 0
    return {"passed": True, "cases": 3, "translation_invariant": True,
            "scope": "Analysis fixtures, without evidence of phase formation or causal support."}


if __name__ == "__main__":
    import json
    print(json.dumps(fixtures()))
