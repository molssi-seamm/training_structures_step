# -*- coding: utf-8 -*-

"""A refused Hessian's replacement is checked too, and a replacement that fails
refuses the molecule (science 2026-11-26: methoxy's analytic UKS Hessian was
refused in SEAMM job 6090 and its finite-difference replacement never checked)."""

import contextlib
import types

import numpy as np
import pytest

import training_structures_step

PASS = {"errors": [0.001, 0.002, 0.001], "max frequency": 3000.0,
        "refused": False, "reason": None, "reference guess": True}  # fmt: skip
FAIL = {"errors": [0.03, 0.01, 0.08], "max frequency": 2989.0, "refused": True,
        "reason": "it disagrees with the gradients by 0.080 (> 0.02)",
        "reference guess": True}  # fmt: skip


class Session:
    def __init__(self, checks):
        self.checks = list(checks)
        self.fd_info = None
        self.checked = []

    def hessian(self, finite_differences=False):
        if finite_differences:
            self.fd_info = {"reference guess": True, "tight SCF": True}
            return np.eye(6), "finite differences over MDI"
        return 2 * np.eye(6), "analytic"

    def check_hessian(self, H, rng, coords_ang=None, **limits):
        self.checked.append((H[0, 0], limits))
        return self.checks.pop(0)


def driver_with(session):
    @contextlib.contextmanager
    def open_session(configuration):
        yield session

    return types.SimpleNamespace(session=open_session)


def run(checks):
    node = training_structures_step.TrainingStructures()
    P = node.parameters.current_values_to_dict(context={})  # quantities, as run()
    P["check the hessian"] = "always"
    session = Session(checks)
    configuration = types.SimpleNamespace(spin_multiplicity=2)
    x0 = np.zeros((2, 3))
    H, info = node._hessian(
        P, driver_with(session), configuration, x0, np.ones(2), np.random.default_rng(0)
    )
    return H, info, session


def test_a_sound_hessian_is_checked_once():
    H, info, session = run([PASS])
    assert info["route"] == "analytic" and H[0, 0] == 2.0
    assert len(session.checked) == 1 and "replacement check" not in info
    assert info["check from reference guess"] is True


def test_a_refused_hessian_is_replaced_and_the_replacement_checked():
    H, info, session = run([FAIL, PASS])
    assert H[0, 0] == 1.0 and info["route"] == "finite differences over MDI"
    assert info["refused"].startswith("it disagrees")
    assert info["replacement check"]["reason"] is None
    assert info["finite differences"]["reference guess"] is True
    # both checks used the step's limits
    assert [h for h, _ in session.checked] == [2.0, 1.0]
    assert session.checked[1][1]["tolerance"] == pytest.approx(0.02)


def test_a_refused_replacement_refuses_the_molecule():
    with pytest.raises(RuntimeError, match="refused too.*not used"):
        run([FAIL, FAIL])
