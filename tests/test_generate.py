# -*- coding: utf-8 -*-

"""The geometry generators and filters (generate.py)."""

import math

import numpy as np
import pytest

from training_structures_step import generate as g

pytest.importorskip("rdkit")


def molecule(smiles, seed=7):
    """symbols, xyz, bonds, orders, conjugated of an MMFF structure."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    Chem.Kekulize(m, clearAromaticFlags=True)
    symbols = [a.GetSymbol() for a in m.GetAtoms()]
    xyz = m.GetConformer().GetPositions()
    bonds, orders, conjugated = [], {}, set()
    for b in m.GetBonds():
        key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        bonds.append(key)
        orders[key] = int(b.GetBondTypeAsDouble())
        if b.GetIsConjugated():
            conjugated.add(key)
    return symbols, xyz, bonds, orders, conjugated


# ----------------------------------------------------------------- connectivity
def test_connectivity_of_a_molecule_matches_its_bonds():
    s, x, bonds, _, _ = molecule("CCO")
    assert g.connectivity(s, x) == set(bonds)
    assert g.connectivity_change(bonds, s, x) is None


def test_a_broken_bond_is_reported():
    s, x, bonds, _, _ = molecule("CO")
    h = next(j for i, j in bonds if s[i] == "O" and s[j] == "H")
    o = next(i for i, j in bonds if j == h)
    y = x.copy()
    y[h] += 1.5 * (x[h] - x[o]) / np.linalg.norm(x[h] - x[o])
    change = g.connectivity_change(bonds, s, y)
    assert change == f"-O{o + 1}-H{h + 1}"


# ----------------------------------------------------------------------- rotors
@pytest.mark.parametrize(
    "smiles, n",
    [
        ("CC", 0),  # methyl rotations are not heavy rotors
        ("CCO", 0),  # nor an OH
        ("CCCC", 1),  # the central C-C
        ("CCCCC", 2),
        ("C1CCCCC1", 0),  # ring bonds are not rotors
        ("CC(=O)NC", 0),  # the amide C-N is conjugated; the others end in methyls
        ("CC(=O)NCC", 1),  # N-ethylacetamide: the N-CH2 bond
        ("CC(=O)OCC", 1),  # ester: only O-CH2 (C(=O)-O is conjugated)
        ("c1ccccc1C(=O)C", 0),  # aryl-C(O) is conjugated
    ],
)
def test_heavy_rotors(smiles, n):
    s, _, bonds, orders, conjugated = molecule(smiles)
    assert len(g.find_rotors(s, bonds, orders, conjugated)) == n


def test_rotors_rank_by_heavy_atoms_moved():
    s, _, bonds, orders, conjugated = molecule("CCCCCC")
    rotors = g.find_rotors(s, bonds, orders, conjugated)
    assert [r["heavy moved"] for r in rotors] == [3, 2, 2]
    assert len(g.find_rotors(s, bonds, orders, conjugated, max_rotors=2)) == 2


def test_set_dihedral_turns_only_the_moving_side():
    s, x, bonds, orders, conjugated = molecule("CCCC")
    (rotor,) = g.find_rotors(s, bonds, orders, conjugated)
    a, i, j, b = rotor["dihedral"]
    for target in (-150.0, -60.0, 0.0, 75.0, 180.0):
        y = g.set_dihedral(x, rotor, target)
        assert g.wrap(g.dihedral(y, a, i, j, b) - target) < 1e-6
        fixed = [k for k in range(len(s)) if k not in rotor["moving"]]
        assert np.allclose(y[fixed], x[fixed])
        # bond lengths are unchanged
        for p, q in bonds:
            assert np.linalg.norm(y[p] - y[q]) == pytest.approx(
                np.linalg.norm(x[p] - x[q])
            )


def test_scan_angles():
    angles = g.scan_angles(60.0, 15.0)
    assert len(angles) == 23 and 60.0 not in angles
    assert angles[0] == 75.0 and angles[-1] == 45.0
    assert all(-180.0 <= a < 180.0 for a in angles)


# ------------------------------------------------------------- bond distortions
def test_bond_classes_take_the_median_bond():
    s, x, bonds, _, _ = molecule("CCO")
    classes = g.bond_classes(s, x, bonds)
    assert set(classes) == {"C-C", "C-H", "C-O", "H-O"}
    i, j = classes["C-H"]
    lengths = sorted(
        np.linalg.norm(x[p] - x[q]) for p, q in bonds if {s[p], s[q]} == {"C", "H"}
    )
    assert np.linalg.norm(x[i] - x[j]) == pytest.approx(
        lengths[(len(lengths) - 1) // 2]
    )
    assert set(g.bond_classes(s, x, bonds, pairs={"H-C", "O-C"})) == {"C-H", "C-O"}


def diatomic(k=0.37, r=0.7414):
    """H2 with a harmonic stretch constant k (E_h/bohr²)."""
    x = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, r]])
    v = g.stretch_vector(x, 0, 1)
    H = k * np.outer(v, v) / (v @ v) ** 2  # vᵀHv = k
    return x, H, k


def test_stretch_amplitude_is_exact_for_a_harmonic_bond():
    x, H, k = diatomic()
    e = 8.0 / g.KCAL_HARTREE
    a = g.stretch_amplitude(H, x, 0, 1, e)
    assert 0.5 * k * (a / g.BOHR) ** 2 == pytest.approx(e)
    assert g.stretch_amplitude(H, x, 0, 1, e, maximum=0.05) == 0.05
    y = g.stretch(x, 0, 1, a)
    assert np.linalg.norm(y[1] - y[0]) == pytest.approx(0.7414 + a)
    assert np.allclose(0.5 * (y[0] + y[1]), 0.5 * (x[0] + x[1]))


def test_diatomic_grid():
    x, H, k = diatomic()
    e = 9.5 / g.KCAL_HARTREE
    grid = g.diatomic_grid(H, x, e, 100)
    assert len(grid) == 100
    a = g.stretch_amplitude(H, x, 0, 1, e)
    assert grid[0] == pytest.approx(-a) and grid[-1] == pytest.approx(a)
    assert len(g.diatomic_grid(H, x, e, 11)) == 10  # r_e itself is left out


# ---------------------------------------------------------------- Hessian check
def harmonic(n_atoms=3, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(3 * n_atoms, 3 * n_atoms))
    K = 0.05 * (A @ A.T) / n_atoms + 0.02 * np.eye(3 * n_atoms)  # E_h/bohr²
    x0 = rng.normal(size=(n_atoms, 3))

    def gradient(x):
        return K @ ((np.asarray(x).reshape(-1) - x0.reshape(-1)) / g.BOHR)

    return K, x0, gradient


def test_a_sound_hessian_passes():
    K, x0, gradient = harmonic()
    r = g.check_hessian(K, gradient, x0, [12.0] * 3, np.random.default_rng(1))
    assert not r["refused"] and max(r["errors"]) < 1e-10 and len(r["errors"]) == 3


def test_a_wrong_hessian_is_refused():
    K, x0, gradient = harmonic()
    bad = K.copy()
    bad[0, 0] *= 1.5
    r = g.check_hessian(bad, gradient, x0, [12.0] * 3, np.random.default_rng(1))
    assert r["refused"] and "disagrees" in r["reason"]


def test_an_absurd_frequency_is_refused():
    K, x0, gradient = harmonic()
    bad = K.copy()
    bad[0, 0] = -1.0e3
    r = g.check_hessian(bad, gradient, x0, [12.0] * 3, np.random.default_rng(1))
    assert r["refused"] and "frequency" in r["reason"]


def test_frequency_units():
    # H2: k = 0.37 E_h/bohr², reduced mass 0.504 amu -> about 4400 cm⁻¹
    x, H, k = diatomic()
    nu = g.frequencies(H, [1.008, 1.008])
    assert max(nu) == pytest.approx(5140.48767 * math.sqrt(k / 0.504), rel=1e-6)


def test_a_stretched_bond_is_exempt_from_the_connectivity_check():
    s, x = ["H", "H"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    y = g.stretch(x, 0, 1, 0.2)  # 0.94 Å: beyond 1.3 x (0.31 + 0.31)
    assert g.connectivity_change([(0, 1)], s, y) == "-H1-H2"
    assert g.connectivity_change([(0, 1)], s, y, exempt=(1, 0)) is None
