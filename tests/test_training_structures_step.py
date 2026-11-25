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
