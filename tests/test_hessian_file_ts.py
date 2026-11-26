# -*- coding: utf-8 -*-

"""Training Structures keeps the Hessian (Paul 2026-11-26): hessian.npz in the
step's directory, its path in the "hessian" record; the export copies it next to
the molecule's .extxyz."""

import types

import numpy as np

from seamm_exec.hessian import load_hessian
from training_structures_step.export_training_set import ExportTrainingSet
from training_structures_step.training_structures import TrainingStructures


def test_the_hessian_is_saved_and_copied(tmp_path):
    root = tmp_path / "job"
    step = root / "5" / "water.c0" / "1"
    step.mkdir(parents=True)
    H = np.eye(9) * 0.3
    x0 = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])
    me = types.SimpleNamespace(flowchart=types.SimpleNamespace(root_directory=root))
    driver = types.SimpleNamespace(
        model_chemistry=lambda: {"level": "ORCA:CMP@R2SCAN-3C"}
    )
    reference = types.SimpleNamespace(charge=0, spin_multiplicity=1)
    info = {"route": "analytic over MDI", "settings": "sha256:6b46", "checked": False}
    name = TrainingStructures._save_hessian(
        me, step, H, [16.0, 1.0, 1.0], x0, ["O", "H", "H"], reference, info, driver
    )
    assert name == "5/water.c0/1/hessian.npz"
    data = load_hessian(root / name)
    assert np.array_equal(data["hessian"], H)
    assert data["info"]["route"] == "analytic over MDI"
    assert data["info"]["model"] == "ORCA:CMP@R2SCAN-3C"
    assert data["info"]["settings"] == "sha256:6b46"
    assert data["info"]["multiplicity"] == 1

    out = root / "training_set"
    out.mkdir()
    copied = ExportTrainingSet._copy_hessian(me, {"file": name}, out, "water")
    assert copied == "water.hessian.npz"
    assert np.array_equal(load_hessian(out / copied)["hessian"], H)
    # an older record, without a file: nothing to copy
    assert ExportTrainingSet._copy_hessian(me, {"route": "x"}, out, "w") is None
