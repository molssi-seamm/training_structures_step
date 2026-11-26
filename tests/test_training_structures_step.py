#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `training_structures_step` package."""

import pytest  # noqa: F401
import training_structures_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = training_structures_step.TrainingStructures()
    assert (
        str(type(result))
        == "<class 'training_structures_step.training_structures.TrainingStructures'>"
    )


def test_hessian_route_default_and_description():
    """The warm engine by default (analytic Hessians, checkable), as the CBH
    approach; separate tasks are described as finite differences."""
    node = training_structures_step.TrainingStructures()
    node._id = (1,)  # a step number for the header, as in a flowchart
    P = node.parameters.values_to_dict()
    assert P["hessian route"] == "warm engine over MDI"
    text = " ".join(node.description_text(P).split())  # unwrapped
    assert "warm engine over MDI" in text and "checked against" in text
    P["hessian route"] = "separate tasks"
    text = " ".join(node.description_text(P).split())
    assert "finite differences of the gradient, run as separate tasks" in text
    assert "checked against" not in text


def test_hessian_route_reaches_the_driver():
    """The chosen route is the HessianDriver's."""
    node = training_structures_step.TrainingStructures()
    for route in ("warm engine over MDI", "automatic", "separate tasks"):
        P = node.parameters.values_to_dict()
        P["hessian route"] = route
        assert node._driver(P).route == route
