# -*- coding: utf-8 -*-

"""A scan point whose torsion ASE cannot hold is dropped, not the molecule
(B1 6091: FixInternals.adjust_positions did not converge, on nitriles)."""

import types

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes

import ase.optimize
import seamm_exec.hessian

from training_structures_step import generate
from training_structures_step.training_structures import TrainingStructures

pytest.importorskip("rdkit")


class Flat(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, engine=None):
        super().__init__()

    def calculate(self, atoms=None, properties=("energy",), changes=all_changes):
        super().calculate(atoms, properties, changes)
        self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}


def test_an_unholdable_torsion_drops_the_point(monkeypatch):
    from .test_generate import molecule

    s, x, bonds, orders, conjugated = molecule("CCCC")
    (rotor,) = generate.find_rotors(s, bonds, orders, conjugated)
    calls = []

    class BFGS:
        def __init__(self, atoms, logfile=None):
            self.atoms = atoms

        def run(self, fmax, steps):
            calls.append(1)
            if len(calls) == 3:
                raise ValueError("FixInternals.adjust_positions did not converge.")
            return True

    monkeypatch.setattr(ase.optimize, "BFGS", BFGS)
    monkeypatch.setattr(seamm_exec.hessian, "_MDICalculator", Flat)
    P = {
        "relaxation convergence": 0.01,
        "relaxation max steps": 10,
        "scan step": 30.0,
    }
    dropped = {}
    points = TrainingStructures._scan(
        types.SimpleNamespace(_as_float=TrainingStructures._as_float),
        P,
        types.SimpleNamespace(engine=None),
        rotor,
        1,
        x,
        s,
        lambda x, energy=None: (None, 0.0),
        dropped,
    )
    assert len(calls) == 11  # 360/30 - 1 points
    assert len(points) == 10
    (why,) = dropped["scan"]
    assert "the torsion could not be held" in why and "did not converge" in why
