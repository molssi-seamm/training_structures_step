# -*- coding: utf-8 -*-

"""Relative energies, the filters, the extended-XYZ file and coverage (export.py)."""

import io
import json

import numpy as np
import pytest

from training_structures_step import export as ex

KCAL = ex.KJ_PER_KCAL
H2O = ["O", "H", "H"]
X0 = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.467], [0.0, -0.757, -0.467]]


def point(name, kind, energy, s2=None, fp="sha256:a", **extra):
    p = {
        "name": name,
        "type": kind,
        "energy": energy,
        "gradients": [[0.1, 0.0, 0.0], [-0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]],
        "fingerprint": fp,
        "S**2": s2,
        "xyz": X0,
    }
    p.update(extra)
    return p


def molecule():
    e0 = -200000.0
    return [
        point("minimum", "minimum", e0),
        point("sample 1", "sample", e0 + 3 * KCAL, **{"target energy": 12.0}),
        point("sample 2", "sample", e0 + 11 * KCAL),  # above a 10 kcal/mol window
        point("scan 1 C2-C3 +15", "scan", e0 + 1.0, **{"scan angle": 15.0}),
        point("stretch O-H +", "stretch", None),  # the label failed
    ]


def test_relative_energies_and_the_window():
    kept, dropped, problem = ex.filter_points(molecule(), window=10 * KCAL)
    assert problem is None
    assert [p["name"] for p in kept] == ["minimum", "sample 1", "scan 1 C2-C3 +15"]
    assert kept[0]["delta energy"] == 0.0
    assert kept[1]["delta energy"] == pytest.approx(3 * KCAL)
    assert dropped == [
        ("sample 2", "sample", "11.00 kcal/mol above the minimum"),
        ("stretch O-H +", "stretch", "no energy (the label failed)"),
    ]
    # half the samples were lost to the window: above the 15 % redraw threshold
    assert ex.sample_drop_fraction(molecule(), dropped) == 0.5


def test_no_window_keeps_everything_labelled():
    kept, dropped, _ = ex.filter_points(molecule())
    assert len(kept) == 4 and len(dropped) == 1


def test_the_spin_check_for_open_shells():
    e0 = -100.0
    points = [
        point("minimum", "minimum", e0, s2=0.7530),
        point("a", "sample", e0 + 1, s2=0.7560),
        point("b", "sample", e0 + 1, s2=0.8200),  # spin contaminated
        point("c", "sample", e0 + 1, s2=None),  # no <S**2>: nothing to compare
    ]
    kept, dropped, _ = ex.filter_points(points, spin_tolerance=0.05)
    assert [p["name"] for p in kept] == ["minimum", "a", "c"]
    assert dropped[0][0] == "b" and "<S**2> 0.8200" in dropped[0][2]


@pytest.mark.parametrize(
    "points, message",
    [
        ([point("a", "sample", -1.0)], "0 points are marked as the minimum"),
        ([point("m", "minimum", None)], "the minimum has no energy"),
        (
            [point("m", "minimum", -1.0), point("a", "sample", -0.9, fp="sha256:b")],
            "2 different settings",
        ),
    ],
)
def test_a_molecule_that_cannot_be_exported(points, message):
    kept, dropped, problem = ex.filter_points(points)
    assert kept == [] and message in problem


def test_the_file_reads_back_with_ase():
    ase_io = pytest.importorskip("ase.io")
    kept, _, _ = ex.filter_points(molecule(), window=10 * KCAL)
    text = ex.molecule_file(H2O, kept, "water = test", units="eV")
    frames = ase_io.read(io.StringIO(text), index=":", format="extxyz")
    assert len(frames) == 3
    first = frames[0]
    # REF_ keys, not a calculator's results; no charge column
    assert first.calc is None
    assert first.info["REF_energy"] == pytest.approx(-200000.0 / ex.KJ_PER_EV)
    assert first.info["delta_energy"] == 0.0
    assert first.info["molecule"] == "water = test"  # quoted, '=' kept
    assert first.info["config_type"] == "minimum"
    forces = first.arrays["REF_forces"]
    assert forces[0][0] == pytest.approx(-0.1 / ex.KJ_PER_EV)  # F = -dE/dx
    assert np.allclose(first.positions, X0)
    assert frames[1].info["target_energy_kJmol"] == 12.0
    assert frames[2].info["scan_angle"] == 15.0
    assert frames[2].info["point"] == "scan 1 C2-C3 +15"


def test_full_precision_positions():
    p = point("minimum", "minimum", -1.0)
    p["delta energy"] = 0.0
    p["xyz"] = [[0.123456789012, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    text = ex.molecule_file(H2O, [p], "w")
    assert "0.123456789012" in text


def test_provenance_record():
    kept, dropped, _ = ex.filter_points(molecule(), window=10 * KCAL)
    record = ex.provenance(
        "water",
        "ORCA:DFT@r2SCAN-3c",
        "sha256:a",
        {"program": "ORCA 6.1.1"},
        {"route": "analytic"},
        42,
        kept,
        dropped,
    )
    assert record["kept"] == {"minimum": 1, "sample": 1, "scan": 1}
    assert record["dropped"][0]["point"] == "sample 2"
    assert json.loads(ex.dumps(record))["random seed"] == 42


def test_coverage_counts_elements_and_bonded_pairs():
    coverage = {}
    xyz = [X0, [[0.0, 0.0, 0.117], [0.0, 0.857, -0.467], [0.0, -0.757, -0.467]]]
    ex.add_coverage(coverage, H2O, [(0, 1), (0, 2)], xyz)
    ex.add_coverage(coverage, ["O", "H", "H"], [(0, 1), (0, 2)], xyz[:1])
    assert coverage["elements"]["O"] == [2, 3] and coverage["elements"]["H"] == [2, 3]
    m, n, r0, r1 = coverage["pairs"]["H-O"]
    assert (m, n) == (2, 3)
    assert r0 == pytest.approx(np.linalg.norm(np.subtract(X0[1], X0[0])))
    assert r1 > r0
    rows = ex.coverage_table(coverage, element_minimum=(2, 3), pair_minimum=(3, 1))
    assert ("element", "H", 2, 3, None, None, True) in rows
    assert [r for r in rows if r[0] == "pair"][0][-1] is False  # 2 < 3 molecules
