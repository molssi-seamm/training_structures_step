# -*- coding: utf-8 -*-

"""The Export Training Set node against a real molsystem database: finding the
training systems, reading the labels, and marking the points."""

import types

import pytest

molsystem = pytest.importorskip("molsystem")

from training_structures_step import ExportTrainingSet  # noqa: E402
from training_structures_step import export as ex  # noqa: E402

MODEL = "ORCA:DFT@r2SCAN-3c"
X0 = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.467], [0.0, -0.757, -0.467]]


def put(conf, name, value, _type="float", units=None):
    if not conf.properties.exists(name):
        conf.properties.add(name, _type, units=units, noerror=True)
    conf.properties.put(name, value)


@pytest.fixture()
def db(tmp_path):
    db = molsystem.SystemDB(filename=str(tmp_path / "seamm.db"))
    other = db.create_system(name="water")
    other.create_configuration().atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol=["O"])
    system = db.create_system(name="water training")
    for k, (name, kind, energy) in enumerate(
        [
            ("minimum", "minimum", -200000.0),
            ("sample 1", "sample", -199990.0),
            ("sample 2", "sample", -199950.0),  # 50 kJ/mol: above the window
        ]
    ):
        conf = system.create_configuration(name=name)
        conf.atoms.append(
            x=[p[0] for p in X0],
            y=[p[1] for p in X0],
            z=[p[2] + 0.01 * k for p in X0],
            symbol=["O", "H", "H"],
        )
        put(conf, "point type#TrainingStructures#generation", kind, "str")
        put(conf, f"energy#Energy#{MODEL}", energy, "float", "kJ/mol")
        put(conf, f"gradients#Energy#{MODEL}", [[0.1, 0, 0]] * 3, "json")
        put(conf, f"settings fingerprint#Energy#{MODEL}", "sha256:a", "str")
        put(conf, f"S**2#Energy#{MODEL}", 0.0, "float")
        if kind == "minimum":
            put(conf, "random seed#TrainingStructures#generation", 42, "int")
            put(conf, "hessian#TrainingStructures#generation", {"route": "x"}, "json")
    return db


def test_points_and_systems(db):
    node = types.SimpleNamespace(
        _is_training_system=ExportTrainingSet._is_training_system
    )
    P = {"systems": "all training systems", "model chemistry": "from the labels"}
    systems = ExportTrainingSet._systems(node, P, db)
    assert [s.name for s in systems] == ["water training"]
    P["systems"] = "water*"
    assert len(ExportTrainingSet._systems(node, P, db)) == 2

    points, model, extra = ExportTrainingSet._points(node, P, systems[0])
    assert model == MODEL
    assert [p["type"] for p in points] == ["minimum", "sample", "sample"]
    assert points[1]["energy"] == -199990.0
    assert points[0]["gradients"] == [[0.1, 0, 0]] * 3
    assert points[0]["S**2"] == 0.0  # the glob characters in the name are exact
    assert points[2]["xyz"][0][2] == pytest.approx(0.137)
    assert extra["seed"] == 42 and extra["hessian"] == {"route": "x"}
    assert extra["symbols"] == ["O", "H", "H"] and extra["multiplicity"] == 1

    kept, dropped, problem = ex.filter_points(points, window=10 * ex.KJ_PER_KCAL)
    assert problem is None and len(kept) == 2 and dropped[0][0] == "sample 2"


def test_points_are_marked_not_deleted(db):
    system = [s for s in db.systems if s.name == "water training"][0]
    conf = system.configurations[2]
    ExportTrainingSet._mark(None, conf, "dropped: too high", None)
    ExportTrainingSet._mark(None, system.configurations[1], "kept", 10.0)
    from training_structures_step.export_training_set import _get

    assert _get(conf, "training status#TrainingStructures#export") == (
        "dropped: too high"
    )
    assert (
        _get(system.configurations[1], "delta energy#TrainingStructures#export") == 10.0
    )
    assert system.n_configurations == 3


def test_file_stem():
    assert ExportTrainingSet._file_stem("C1=CC=C1 training") == "C1_CC_C1_training"
