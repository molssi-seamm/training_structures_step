# -*- coding: utf-8 -*-

"""The labelled training points of one molecule: relative energies, the filters,
the extended-XYZ file and the coverage counts, as plain functions (no flowchart
state) so that they can be tested on their own.

Energies are in kJ/mol and gradients in kJ/mol/Å, as the Energy step stores them,
unless a name says otherwise.
"""

import json
import math

KJ_PER_EV = 96.48533212331002
KJ_PER_HARTREE = 2625.4996394799
KJ_PER_KCAL = 4.184

#: the energy and force units of the file: name -> (kJ/mol per unit, force units)
UNITS = {
    "eV": (KJ_PER_EV, "eV/Å"),
    "kJ/mol": (1.0, "kJ/mol/Å"),
    "E_h": (KJ_PER_HARTREE, "E_h/Å"),
}


def filter_points(
    points,
    window=None,
    spin_tolerance=None,
):
    """Relative energies and the filters, for one molecule's labelled points.

    Parameters
    ----------
    points : [dict]
        One per configuration, in order, each with "name", "type" (the Training
        Structures point type), "energy" (kJ/mol or None when unlabelled),
        "gradients", "fingerprint" (the settings fingerprint or None) and "S**2"
        (or None). Exactly one has type "minimum".
    window : float or None
        Drop a point more than this (kJ/mol) above the minimum.
    spin_tolerance : float or None
        For an open shell, drop a point whose <S**2> differs from the minimum's by
        more than this.

    Returns
    -------
    (kept, dropped, problem)
        kept: the points with "delta energy" added, the minimum first; dropped:
        (name, type, reason); problem: None, or why the molecule as a whole cannot
        be exported (no labelled minimum, mixed settings).
    """
    minima = [p for p in points if p.get("type") == "minimum"]
    if len(minima) != 1:
        return [], [], f"{len(minima)} points are marked as the minimum"
    minimum = minima[0]
    if minimum.get("energy") is None:
        return [], [], "the minimum has no energy"
    fingerprints = {p["fingerprint"] for p in points if p.get("energy") is not None}
    fingerprints.discard(None)
    if len(fingerprints) > 1:
        return (
            [],
            [],
            (
                f"the points were labelled with {len(fingerprints)} different "
                "settings; one molecule must not mix levels or settings"
            ),
        )
    e0 = float(minimum["energy"])
    s0 = minimum.get("S**2")
    kept, dropped = [], []
    for p in [minimum] + [q for q in points if q is not minimum]:
        if p.get("energy") is None:
            dropped.append((p["name"], p.get("type"), "no energy (the label failed)"))
            continue
        de = float(p["energy"]) - e0
        if window is not None and de > window:
            dropped.append(
                (
                    p["name"],
                    p.get("type"),
                    f"{de / KJ_PER_KCAL:.2f} kcal/mol above the minimum",
                )
            )
            continue
        s2 = p.get("S**2")
        if spin_tolerance is not None and s0 is not None and s2 is not None:
            if abs(s2 - s0) > spin_tolerance:
                dropped.append(
                    (
                        p["name"],
                        p.get("type"),
                        f"<S**2> {s2:.4f} against {s0:.4f} at the minimum",
                    )
                )
                continue
        kept.append(dict(p, **{"delta energy": de}))
    return kept, dropped, None


def merge_redraw(points, redraw_points, same_model=True, tolerance=0.01):
    """A molecule's points with a redraw's samples in place of the first pass's.

    The first pass's samples are kept in the list, marked "replaced" (so the
    caller can mark their configurations); the redraw's minimum is a duplicate of
    the first pass's and is left out. Both minima must agree in energy (within
    ``tolerance`` kJ/mol) and settings: the redraw must start from the same
    structure at the same level.

    Returns
    -------
    (points, problem)
    """
    if not same_model:
        return points, "it was labelled with another model chemistry"
    first_min = next((p for p in points if p.get("type") == "minimum"), None)
    redraw_min = next((p for p in redraw_points if p.get("type") == "minimum"), None)
    if first_min is None or redraw_min is None:
        return points, "a pass has no minimum"
    if first_min.get("energy") is None or redraw_min.get("energy") is None:
        return points, "a minimum has no energy"
    if first_min.get("fingerprint") != redraw_min.get("fingerprint"):
        return points, "its settings differ from the first pass's"
    if abs(float(first_min["energy"]) - float(redraw_min["energy"])) > tolerance:
        return points, "it started from a different minimum"
    merged = []
    for p in points:
        if p.get("type") == "sample":
            merged.append(dict(p, replaced=True))
        else:
            merged.append(p)
    for p in redraw_points:
        if p.get("type") == "sample":
            merged.append(dict(p, name=f"redraw {p['name']}"))
    return merged, None


def sample_drop_fraction(points, dropped):
    """The fraction of the normal-mode samples dropped by the window or the spin
    check (the recipe redraws a molecule's samples when it exceeds 15 %)."""
    n = sum(1 for p in points if p.get("type") == "sample")
    if n == 0:
        return 0.0
    lost = sum(
        1
        for name, kind, reason in dropped
        if kind == "sample" and not reason.startswith("no energy")
    )
    return lost / n


def _value(value):
    """An extended-XYZ key value: numbers as they are, text quoted when needed."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    text = str(value)
    if text == "" or any(c in text for c in ' ="\t,'):
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return text


def _key(name):
    """A key without spaces or '=' (ASE would split or merge such keys)."""
    return "".join(c if (c.isalnum() or c == "_") else "_" for c in str(name))


def extxyz_frame(symbols, xyz, info, gradients=None, units="eV"):
    """One extended-XYZ frame: species, positions (full precision) and, with
    gradients (kJ/mol/Å), the forces in the file's units as REF_forces.

    ``info`` gives the comment-line keys; numbers already in the file's units.
    The energy goes in as REF_energy (not 'energy', which ASE would turn into a
    calculator result), and there is no 'charge' column (ASE reads that as a
    calculator result too)."""
    per_kj = 1.0 / UNITS[units][0]
    properties = "species:S:1:pos:R:3" + (":REF_forces:R:3" if gradients else "")
    items = [f"Properties={properties}", 'pbc="F F F"']
    items += [f"{_key(k)}={_value(v)}" for k, v in info.items() if v is not None]
    lines = [str(len(symbols)), " ".join(items)]
    for i, (s, (x, y, z)) in enumerate(zip(symbols, xyz)):
        row = f"{s:<2s} {x:20.12f} {y:20.12f} {z:20.12f}"
        if gradients:
            fx, fy, fz = (-g * per_kj for g in gradients[i])
            row += f" {fx:20.12f} {fy:20.12f} {fz:20.12f}"
        lines.append(row)
    return "\n".join(lines) + "\n"


#: point keys copied into each frame's comment line, with their file key
FRAME_KEYS = (
    ("type", "config_type"),
    ("name", "point"),
    ("target energy", "target_energy_kJmol"),
    ("harmonic energy", "harmonic_energy_kJmol"),
    ("prescreen energy", "prescreen_energy_kJmol"),
    ("scan angle", "scan_angle"),
    ("bond class", "bond_class"),
    ("bond change", "bond_change"),
    ("S**2", "S2"),
)


def molecule_file(symbols, kept, molecule, units="eV", charge=0, multiplicity=1):
    """The extended-XYZ text of one molecule's kept points (the minimum first)."""
    per_kj = 1.0 / UNITS[units][0]
    frames = []
    for p in kept:
        info = {
            "REF_energy": p["energy"] * per_kj,
            "delta_energy": p["delta energy"] * per_kj,
            "energy_units": units,
            "molecule": molecule,
            "total_charge": int(charge),
            "multiplicity": int(multiplicity),
        }
        for key, out in FRAME_KEYS:
            if p.get(key) is not None:
                info[out] = p[key]
        frames.append(extxyz_frame(symbols, p["xyz"], info, p.get("gradients"), units))
    return "".join(frames)


def provenance(molecule, model, fingerprint, record, hessian, seed, kept, dropped):
    """The sidecar record of one molecule's export (JSON-ready)."""
    counts = {}
    for p in kept:
        counts[p["type"]] = counts.get(p["type"], 0) + 1
    return {
        "molecule": molecule,
        "model chemistry": model,
        "settings fingerprint": fingerprint,
        "calculation": record,
        "hessian": hessian,
        "random seed": seed,
        "kept": counts,
        "dropped": [
            {"point": name, "type": kind, "reason": reason}
            for name, kind, reason in dropped
        ],
    }


def dumps(record):
    return json.dumps(record, indent=1, sort_keys=False, default=str)


# --------------------------------------------------------------------------- #
# Coverage
# --------------------------------------------------------------------------- #
def add_coverage(coverage, symbols, bonds, xyz_list):
    """Count one molecule's kept points: every point counts for every element and
    every bonded element pair of the molecule (the E2a census rule), and each
    pair's sampled bond lengths are tracked.

    coverage : {"elements": {el: [molecules, points]}, "pairs": {pair: [molecules,
    points, r_min, r_max]}}, updated in place."""
    n = len(xyz_list)
    if n == 0:
        return coverage
    elements = coverage.setdefault("elements", {})
    pairs = coverage.setdefault("pairs", {})
    for el in sorted(set(symbols)):
        entry = elements.setdefault(el, [0, 0])
        entry[0] += 1
        entry[1] += n
    seen = set()
    for i, j in bonds:
        key = "-".join(sorted((symbols[i], symbols[j])))
        entry = pairs.setdefault(key, [0, 0, math.inf, -math.inf])
        if key not in seen:
            seen.add(key)
            entry[0] += 1
            entry[1] += n
        for xyz in xyz_list:
            r = math.dist(xyz[i], xyz[j])
            entry[2] = min(entry[2], r)
            entry[3] = max(entry[3], r)
    return coverage


def coverage_table(coverage, element_minimum=(30, 500), pair_minimum=(10, 100)):
    """Rows (kind, name, molecules, points, r_min, r_max, passes) of a coverage
    record, against the minimum (molecules, points) for elements and pairs."""
    rows = []
    for el, (m, n) in sorted(coverage.get("elements", {}).items()):
        ok = m >= element_minimum[0] and n >= element_minimum[1]
        rows.append(("element", el, m, n, None, None, ok))
    for pair, (m, n, r0, r1) in sorted(coverage.get("pairs", {}).items()):
        ok = m >= pair_minimum[0] and n >= pair_minimum[1]
        rows.append(("pair", pair, m, n, r0, r1, ok))
    return rows
