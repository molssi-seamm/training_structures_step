# -*- coding: utf-8 -*-

"""The geometry generators and filters of the Training Structures step, as plain
functions of coordinates, bonds and Hessians (no flowchart state), so that they
can be tested on their own.

Units: coordinates in Å, Hessians in E_h/bohr², gradients in E_h/bohr, energies in
E_h unless a name says otherwise.
"""

import math

import numpy as np

BOHR = 0.529177210903  # Å
KCAL_HARTREE = 627.509474  # kcal/mol per E_h
KJ_HARTREE = 2625.4996394799  # kJ/mol per E_h


# --------------------------------------------------------------------------- #
# Connectivity
# --------------------------------------------------------------------------- #
def covalent_radii(symbols):
    """RDKit's covalent radii (Å) of the elements."""
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    return np.array([table.GetRcovalent(s) for s in symbols])


def connectivity(symbols, xyz, factor=1.3, radii=None):
    """The bonds (i < j) of a structure: atoms closer than ``factor`` times the
    sum of their covalent radii."""
    x = np.asarray(xyz, dtype=float).reshape(-1, 3)
    r = covalent_radii(symbols) if radii is None else np.asarray(radii)
    bonds = set()
    for i in range(len(x)):
        d = np.linalg.norm(x[i + 1 :] - x[i], axis=1)
        for k in np.nonzero(d < factor * (r[i + 1 :] + r[i]))[0]:
            bonds.add((i, i + 1 + int(k)))
    return bonds


def connectivity_change(
    reference_bonds, symbols, xyz, factor=1.3, radii=None, exempt=None
):
    """None if the structure has exactly the reference bonds, else the bonds
    broken and formed, as text (e.g. '-C1-H5 +H5-O3'). ``exempt`` is a bond
    (i, j) not judged: a bond stretched on purpose."""
    now = connectivity(symbols, xyz, factor, radii)
    ref = {tuple(sorted(b)) for b in reference_bonds}
    if exempt is not None:
        key = tuple(sorted(exempt))
        now.discard(key)
        ref.discard(key)
    if now == ref:
        return None

    def label(b):
        i, j = b
        return f"{symbols[i]}{i + 1}-{symbols[j]}{j + 1}"

    return " ".join(
        [f"-{label(b)}" for b in sorted(ref - now)]
        + [f"+{label(b)}" for b in sorted(now - ref)]
    )


# --------------------------------------------------------------------------- #
# Rotors and torsion scans
# --------------------------------------------------------------------------- #
def _side(n_atoms, bonds, i, j):
    """The atoms on j's side of the bond i-j (j included), or None when the bond
    is in a ring."""
    neighbours = [set() for _ in range(n_atoms)]
    for a, b in bonds:
        neighbours[a].add(b)
        neighbours[b].add(a)
    seen, todo = {j}, [j]
    while todo:
        a = todo.pop()
        for b in neighbours[a]:
            if a == j and b == i:
                continue
            if b == i:
                return None  # back to i another way: a ring
            if b not in seen:
                seen.add(b)
                todo.append(b)
    return seen


def find_rotors(symbols, bonds, orders=None, conjugated=None, max_rotors=None):
    """The heavy rotors of a molecule, most heavy atoms moved first.

    A heavy rotor is an acyclic single bond i-j, not conjugated, where i and j
    each have another heavy-atom neighbour (so methyl, OH and NH2 rotations are
    not rotors). Rotating about it moves the smaller side; rotors are ranked by
    the number of heavy atoms that moves, ties going to the lower atom indices.

    Parameters
    ----------
    symbols : [str]
    bonds : [(int, int)]
        0-based.
    orders : {(i, j): int} or None
        Bond orders (i < j); a bond without one counts as single.
    conjugated : set of (i, j) or None
        The conjugated bonds (i < j), e.g. from RDKit's GetIsConjugated.
    max_rotors : int or None

    Returns
    -------
    [dict]
        Each with "bond" (i, j), "dihedral" (a, i, j, b), "moving" (the atoms
        that rotate, a sorted list) and "heavy moved".
    """
    n = len(symbols)
    heavy = [s != "H" for s in symbols]
    neighbours = [set() for _ in range(n)]
    for a, b in bonds:
        neighbours[a].add(b)
        neighbours[b].add(a)
    orders = {} if orders is None else orders
    conjugated = set() if conjugated is None else conjugated
    rotors = []
    for a, b in sorted({tuple(sorted(x)) for x in bonds}):
        if not (heavy[a] and heavy[b]):
            continue
        if orders.get((a, b), 1) != 1 or (a, b) in conjugated:
            continue
        ha = sorted(x for x in neighbours[a] - {b} if heavy[x])
        hb = sorted(x for x in neighbours[b] - {a} if heavy[x])
        if not ha or not hb:
            continue
        side_b = _side(n, bonds, a, b)
        if side_b is None:
            continue  # in a ring
        side_a = set(range(n)) - side_b
        heavy_a = sum(1 for x in side_a if heavy[x])
        heavy_b = sum(1 for x in side_b if heavy[x])
        moving = side_b if heavy_b <= heavy_a else side_a
        rotors.append(
            {
                "bond": (a, b),
                "dihedral": (ha[0], a, b, hb[0]),
                "moving": sorted(moving),
                "heavy moved": min(heavy_a, heavy_b),
            }
        )
    rotors.sort(key=lambda r: (-r["heavy moved"], r["bond"]))
    return rotors if max_rotors is None else rotors[:max_rotors]


def dihedral(x, a, b, c, d):
    """The dihedral a-b-c-d (degrees, -180..180)."""
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def rotate_about_bond(xyz, i, j, moving, degrees):
    """A copy of xyz with the atoms ``moving`` turned by ``degrees`` about the
    i-j axis (right-handed about i -> j)."""
    x = np.array(xyz, dtype=float).reshape(-1, 3)
    axis = x[j] - x[i]
    axis /= np.linalg.norm(axis)
    t = math.radians(degrees)
    K = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    R = np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * K @ K
    idx = list(moving)
    x[idx] = (x[idx] - x[j]) @ R.T + x[j]
    return x


def set_dihedral(xyz, rotor, target):
    """A copy of xyz with the rotor's dihedral set to ``target`` degrees by
    turning its moving side."""
    a, i, j, b = rotor["dihedral"]
    x = np.asarray(xyz, dtype=float).reshape(-1, 3)
    now = dihedral(x, a, i, j, b)
    # turning the side that holds b about i -> j by t changes the dihedral by t;
    # turning the side that holds a changes it by -t
    sign = 1.0 if b in rotor["moving"] else -1.0
    return rotate_about_bond(x, i, j, rotor["moving"], sign * (target - now))


def scan_angles(start, step=15.0):
    """The scan's dihedral targets: start + k*step, k = 1 ... 360/step - 1 (the
    start itself is the minimum and is not repeated), wrapped to -180..180."""
    n = int(round(360.0 / step))
    return [((start + k * step + 180.0) % 360.0) - 180.0 for k in range(1, n)]


# --------------------------------------------------------------------------- #
# Bond distortions
# --------------------------------------------------------------------------- #
def bond_classes(symbols, xyz, bonds, pairs=None):
    """{"C-H": (i, j), ...}: for each bonded element pair, the bond of median
    length (the lower of the two middle ones; ties to the lower indices).
    ``pairs`` limits it to those classes (e.g. {"C-H", "C-C"}; order of the
    elements does not matter)."""
    x = np.asarray(xyz, dtype=float).reshape(-1, 3)
    groups = {}
    for i, j in sorted({tuple(sorted(b)) for b in bonds}):
        key = "-".join(sorted((symbols[i], symbols[j])))
        groups.setdefault(key, []).append((float(np.linalg.norm(x[i] - x[j])), i, j))
    if pairs is not None:
        wanted = {"-".join(sorted(p.replace(" ", "").split("-"))) for p in pairs}
        groups = {k: v for k, v in groups.items() if k in wanted}
    out = {}
    for key, items in sorted(groups.items()):
        items.sort()
        _, i, j = items[(len(items) - 1) // 2]
        out[key] = (i, j)
    return out


def stretch_vector(xyz, i, j):
    """The unit-free stretch direction (3N): -u/2 on i, +u/2 on j, u the unit
    vector i -> j; a displacement a times it lengthens the bond by a."""
    x = np.asarray(xyz, dtype=float).reshape(-1, 3)
    u = (x[j] - x[i]) / np.linalg.norm(x[j] - x[i])
    v = np.zeros_like(x)
    v[i] = -0.5 * u
    v[j] = 0.5 * u
    return v.reshape(-1)


def stretch_amplitude(hessian, xyz, i, j, energy, maximum=None):
    """The bond-length change a (Å) whose harmonic energy along the stretch is
    ``energy`` (E_h): 1/2 a² vᵀHv = energy, exact for this two-atom
    displacement; capped at ``maximum`` (Å)."""
    v = stretch_vector(xyz, i, j)
    k = float(v @ np.asarray(hessian, dtype=float) @ v)  # E_h/bohr²
    if k <= 0.0:
        a = maximum
    else:
        a = math.sqrt(2.0 * energy / k) * BOHR
        if maximum is not None:
            a = min(a, maximum)
    return a


def stretch(xyz, i, j, a):
    """A copy of xyz with bond i-j lengthened by a (Å; negative shortens), each
    atom moving a/2."""
    x = np.asarray(xyz, dtype=float).reshape(-1)
    return (x + a * stretch_vector(xyz, i, j)).reshape(-1, 3)


def diatomic_grid(hessian, xyz, energy, n_points):
    """Bond-length changes (Å) for a diatomic's 1-D stretch grid: n_points equally
    spaced from -a to +a, a the displacement whose harmonic energy is ``energy``
    (E_h); r_e itself (the minimum, already a point) is left out if it falls on
    the grid."""
    a = stretch_amplitude(hessian, xyz, 0, 1, energy)
    grid = np.linspace(-a, a, n_points)
    return [float(g) for g in grid if abs(g) > 1e-9 * a]


# --------------------------------------------------------------------------- #
# Hessian check
# --------------------------------------------------------------------------- #
def frequencies(hessian, masses_amu):
    """Signed harmonic frequencies (cm⁻¹) of a Cartesian Hessian (E_h/bohr²),
    without projecting translations and rotations (for the frequency guard)."""
    m = np.repeat(np.asarray(masses_amu, dtype=float), 3)
    H = np.asarray(hessian, dtype=float)
    H_mw = H / np.sqrt(np.outer(m, m))
    lam = np.linalg.eigvalsh(0.5 * (H_mw + H_mw.T))
    # E_h/(bohr² amu) -> (cm⁻¹)²: 5140.487^2 per unit
    factor = 5140.48767
    return np.sign(lam) * np.sqrt(np.abs(lam)) * factor


def check_hessian(
    hessian,
    gradient,
    xyz,
    masses_amu,
    rng,
    n_displacements=3,
    rms=0.005,
    tolerance=0.02,
    max_frequency=5000.0,
):
    """Test a Hessian against central differences of the gradient.

    For each of ``n_displacements`` random displacements δ with an RMS of ``rms``
    Å, e = ‖Hδ − [g(x₀+δ) − g(x₀−δ)]/2‖ / ‖Hδ‖ (2 gradients each). The Hessian
    is refused if any e exceeds ``tolerance`` or any frequency magnitude exceeds
    ``max_frequency`` cm⁻¹.

    Parameters
    ----------
    gradient : callable
        gradient(xyz (N, 3) Å) -> (3N,) E_h/bohr.

    Returns
    -------
    dict
        "errors": [e], "max frequency": the largest frequency magnitude, "refused": bool and
        "reason": text or None.
    """
    x0 = np.asarray(xyz, dtype=float).reshape(-1)
    H = np.asarray(hessian, dtype=float)
    nu = frequencies(H, masses_amu)
    biggest = float(np.max(np.abs(nu)))
    errors = []
    for _ in range(n_displacements):
        d = rng.standard_normal(x0.size)
        d *= rms / math.sqrt(np.mean(d**2))  # Å
        Hd = H @ (d / BOHR)
        g_plus = np.asarray(gradient((x0 + d).reshape(-1, 3)), dtype=float)
        g_minus = np.asarray(gradient((x0 - d).reshape(-1, 3)), dtype=float)
        dg = 0.5 * (g_plus - g_minus)
        errors.append(float(np.linalg.norm(Hd - dg) / np.linalg.norm(Hd)))
    reason = None
    if biggest > max_frequency:
        reason = f"a frequency of {biggest:.0f} cm⁻¹ (> {max_frequency:.0f})"
    elif max(errors) > tolerance:
        reason = f"it disagrees with the gradients by {max(errors):.3f} (> {tolerance})"
    return {
        "errors": errors,
        "max frequency": biggest,
        "refused": reason is not None,
        "reason": reason,
    }


def wrap(delta):
    """The magnitude of an angle difference, wrapped to 0..180 degrees."""
    return abs((delta + 180.0) % 360.0 - 180.0)
